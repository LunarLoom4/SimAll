// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/CanteraParser.cpp
// Phase  : 11
// =============================================================================
#include "combustion/CanteraParser.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace simall::combustion {

namespace {

std::string slurp(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("CanteraParser: cannot open " + path);
    std::ostringstream ss; ss << in.rdbuf();
    return ss.str();
}

CanteraFormat sniff(const std::string& path) {
    auto ext = std::filesystem::path(path).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    if (ext == ".yaml" || ext == ".yml") return CanteraFormat::Yaml;
    if (ext == ".cti")                   return CanteraFormat::Cti;
    return CanteraFormat::Yaml;
}

inline void strip_comment(std::string& s) {
    auto h = s.find('#');
    if (h != std::string::npos) s.erase(h);
}
inline void rstrip(std::string& s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
}
inline std::string trim_copy(std::string s) {
    auto l = std::find_if(s.begin(), s.end(),
                          [](unsigned char c){ return !std::isspace(c); });
    s.erase(s.begin(), l);
    rstrip(s);
    return s;
}

inline std::size_t indent_of(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size() && s[i] == ' ') ++i;
    return i;
}

std::vector<std::string> split_inline(const std::string& v, char open, char close) {
    std::vector<std::string> out;
    auto l = v.find(open);
    auto r = v.rfind(close);
    if (l == std::string::npos || r == std::string::npos || r <= l) return out;
    std::string body = v.substr(l + 1, r - l - 1);
    int depth = 0;
    std::string cur;
    for (char c : body) {
        if      (c == '{' || c == '[') { ++depth; cur += c; }
        else if (c == '}' || c == ']') { --depth; cur += c; }
        else if (c == ',' && depth == 0) { out.push_back(trim_copy(cur)); cur.clear(); }
        else cur += c;
    }
    if (!cur.empty()) out.push_back(trim_copy(cur));
    return out;
}

double parse_double(const std::string& s, double fallback = 0.0) {
    try { return std::stod(s); } catch (...) { return fallback; }
}

void parse_equation(const std::string& eq,
                    std::map<std::string, double>& reactants,
                    std::map<std::string, double>& products,
                    ReactionDirection& dir) {
    std::size_t sep = std::string::npos;
    int seplen = 0;
    if ((sep = eq.find("<=>")) != std::string::npos) {
        seplen = 3; dir = ReactionDirection::Reversible;
    } else if ((sep = eq.find("=>")) != std::string::npos) {
        seplen = 2; dir = ReactionDirection::Forward;
    } else if ((sep = eq.find('=')) != std::string::npos) {
        seplen = 1; dir = ReactionDirection::Reversible;
    } else {
        return;
    }

    auto side = [](const std::string& s, std::map<std::string, double>& out) {
        std::stringstream ss(s);
        std::string tok;
        std::vector<std::string> terms;
        while (std::getline(ss, tok, '+')) terms.push_back(trim_copy(tok));
        for (auto& t : terms) {
            if (t.empty() || t == "M") continue;
            double nu = 1.0;
            std::string species = t;
            std::size_t i = 0;
            while (i < t.size() &&
                   (std::isdigit(static_cast<unsigned char>(t[i])) || t[i] == '.')) ++i;
            if (i > 0 && i < t.size() && std::isspace(static_cast<unsigned char>(t[i]))) {
                nu = parse_double(t.substr(0, i), 1.0);
                species = trim_copy(t.substr(i));
            }
            if (!species.empty()) out[species] += nu;
        }
    };

    side(eq.substr(0, sep),                            reactants);
    side(eq.substr(sep + static_cast<std::size_t>(seplen)), products);
}

}  // namespace

ChemkinMechanism CanteraParser::parse(const std::string& path,
                                       const CanteraReadOptions& opts) {
    auto fmt = (opts.format == CanteraFormat::Auto) ? sniff(path) : opts.format;
    return parseString(slurp(path), fmt, opts);
}

ChemkinMechanism CanteraParser::parseString(const std::string& text,
                                             CanteraFormat fmt,
                                             const CanteraReadOptions& opts) {
    if (fmt == CanteraFormat::Cti) {
        SIMALL_LOG_WARN("combustion.cantera",
            "CTI input: please convert via 'ck2yaml'; treating as YAML.");
    }

    ChemkinMechanism mech;
    mech.energyUnit   = EnergyUnit::JoulesPerMole;
    mech.quantityUnit = QuantityUnit::Moles;

    std::istringstream iss(text);
    std::string raw;

    enum class Section { None, Phases, Species, Reactions, Elements };
    Section section = Section::None;

    ChemkinSpecies  curSpec;   bool inSpecies  = false;
    ChemkinReaction curRx;     bool inReaction = false;

    auto flush_species = [&]() {
        if (inSpecies && !curSpec.name.empty()) mech.species.push_back(std::move(curSpec));
        curSpec = {};
        inSpecies = false;
    };
    auto flush_reaction = [&]() {
        if (inReaction && !curRx.equation.empty()) {
            parse_equation(curRx.equation, curRx.reactants, curRx.products, curRx.direction);
            mech.reactions.push_back(std::move(curRx));
        }
        curRx = {};
        inReaction = false;
    };

    while (std::getline(iss, raw)) {
        std::string line = raw;
        strip_comment(line);
        rstrip(line);
        if (line.empty()) continue;

        const std::string body = trim_copy(line);

        if (indent_of(line) == 0 && body.back() == ':') {
            flush_species();
            flush_reaction();
            const std::string key = body.substr(0, body.size() - 1);
            if      (key == "phases")    section = Section::Phases;
            else if (key == "species")   section = Section::Species;
            else if (key == "reactions") section = Section::Reactions;
            else if (key == "elements")  section = Section::Elements;
            else                          section = Section::None;
            continue;
        }

        if (body.rfind("- ", 0) == 0) {
            if (section == Section::Species)   { flush_species();  inSpecies  = true; }
            if (section == Section::Reactions) { flush_reaction(); inReaction = true; }
            if (section == Section::Elements) {
                ChemkinElement e;
                e.symbol = trim_copy(body.substr(2));
                if (!e.symbol.empty()) mech.elements.push_back(e);
                continue;
            }
            std::string item = body.substr(2);
            auto col = item.find(':');
            if (col != std::string::npos) {
                std::string k = trim_copy(item.substr(0, col));
                std::string v = trim_copy(item.substr(col + 1));
                if (section == Section::Species && k == "name") {
                    curSpec.name = v;
                } else if (section == Section::Reactions && k == "equation") {
                    if (!v.empty() && (v.front() == '"' || v.front() == '\''))
                        v = v.substr(1, v.size() - 2);
                    curRx.equation = v;
                }
            }
            continue;
        }

        if (inSpecies || inReaction) {
            auto col = body.find(':');
            if (col == std::string::npos) continue;
            std::string k = trim_copy(body.substr(0, col));
            std::string v = trim_copy(body.substr(col + 1));

            if (inSpecies) {
                if (k == "composition") {
                    for (const auto& kv : split_inline(v, '{', '}')) {
                        auto c = kv.find(':');
                        if (c == std::string::npos) continue;
                        const std::string el = trim_copy(kv.substr(0, c));
                        const int        n   = static_cast<int>(parse_double(trim_copy(kv.substr(c + 1))));
                        if (!el.empty()) curSpec.composition[el] = n;
                    }
                } else if (k == "thermo") {
                    curSpec.hasThermo   = true;
                    curSpec.thermoIndex = static_cast<int>(mech.species.size());
                }
            }

            if (inReaction) {
                if (k == "rate-constant" || k == "rate-coefficient") {
                    for (const auto& kv : split_inline(v, '{', '}')) {
                        auto c = kv.find(':');
                        if (c == std::string::npos) continue;
                        const std::string kk = trim_copy(kv.substr(0, c));
                        const double      d  = parse_double(trim_copy(kv.substr(c + 1)));
                        if      (kk == "A")  curRx.fwd.A    = d;
                        else if (kk == "b")  curRx.fwd.beta = d;
                        else if (kk == "Ea") curRx.fwd.Ea   = d;
                    }
                } else if (k == "low-P-rate-constant" || k == "low") {
                    curRx.isLindemann = true;
                    for (const auto& kv : split_inline(v, '{', '}')) {
                        auto c = kv.find(':');
                        if (c == std::string::npos) continue;
                        const std::string kk = trim_copy(kv.substr(0, c));
                        const double      d  = parse_double(trim_copy(kv.substr(c + 1)));
                        if      (kk == "A")  curRx.low.A    = d;
                        else if (kk == "b")  curRx.low.beta = d;
                        else if (kk == "Ea") curRx.low.Ea   = d;
                    }
                } else if (k == "Troe") {
                    curRx.troe.enabled = true;
                    for (const auto& kv : split_inline(v, '{', '}')) {
                        auto c = kv.find(':');
                        if (c == std::string::npos) continue;
                        const std::string kk = trim_copy(kv.substr(0, c));
                        const double      d  = parse_double(trim_copy(kv.substr(c + 1)));
                        if      (kk == "A")   curRx.troe.a  = d;
                        else if (kk == "T1")  curRx.troe.T1 = d;
                        else if (kk == "T2")  curRx.troe.T2 = d;
                        else if (kk == "T3")  curRx.troe.T3 = d;
                    }
                } else if (k == "efficiencies") {
                    for (const auto& kv : split_inline(v, '{', '}')) {
                        auto c = kv.find(':');
                        if (c == std::string::npos) continue;
                        const std::string sp = trim_copy(kv.substr(0, c));
                        const double      d  = parse_double(trim_copy(kv.substr(c + 1)));
                        if (!sp.empty()) curRx.thirdBodyEff[sp] = d;
                    }
                } else if (k == "duplicate") {
                    curRx.isDuplicate = (v == "true" || v == "True" || v == "yes");
                }
            }
        }
    }
    flush_species();
    flush_reaction();

    if (opts.stripUnused) {
        std::set<std::string> used;
        for (const auto& r : mech.reactions) {
            for (const auto& kv : r.reactants) used.insert(kv.first);
            for (const auto& kv : r.products)  used.insert(kv.first);
        }
        mech.species.erase(std::remove_if(mech.species.begin(), mech.species.end(),
            [&](const ChemkinSpecies& s){ return used.find(s.name) == used.end(); }),
            mech.species.end());
    }
    if (opts.validateThermo) {
        std::size_t bad = 0;
        for (const auto& s : mech.species) if (!s.hasThermo) ++bad;
        if (bad) {
            SIMALL_LOG_WARN("combustion.cantera",
                bad, " of ", mech.species.size(),
                " species lack NASA-7 thermo blocks (load separately via ThermoNasaParser).");
        }
    }

    SIMALL_LOG_INFO("combustion.cantera",
        "Parsed Cantera mechanism: ", mech.elements.size(), " elements, ",
        mech.species.size(), " species, ", mech.reactions.size(), " reactions.");
    return mech;
}

}  // namespace simall::combustion
