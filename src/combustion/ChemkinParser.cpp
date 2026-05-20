// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/ChemkinParser.cpp
// =============================================================================
#include "combustion/ChemkinParser.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace simall::combustion
{

namespace
{

std::string trim(std::string s)
{
    auto issp = [](unsigned char c) { return std::isspace(c) != 0; };
    while (!s.empty() && issp(s.front()))
        s.erase(s.begin());
    while (!s.empty() && issp(s.back()))
        s.pop_back();
    return s;
}

std::string upper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
    return s;
}

std::vector<std::string> split_ws(const std::string& s)
{
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}

double atomic_weight(const std::string& sym)
{
    // Most-common-isotope average weights (g/mol) for the common elements.
    static const std::map<std::string, double> tbl = {{"H", 1.00794},
                                                      {"D", 2.0141},
                                                      {"C", 12.0107},
                                                      {"N", 14.0067},
                                                      {"O", 15.9994},
                                                      {"AR", 39.948},
                                                      {"HE", 4.0026},
                                                      {"NE", 20.1797},
                                                      {"S", 32.065},
                                                      {"F", 18.9984},
                                                      {"CL", 35.453},
                                                      {"BR", 79.904},
                                                      {"I", 126.9045},
                                                      {"P", 30.9738},
                                                      {"SI", 28.0855}};
    auto it = tbl.find(upper(sym));
    return (it == tbl.end()) ? 0.0 : it->second;
}

// Strip end-of-line comments starting with '!' (Chemkin comment marker).
std::string strip_comment(const std::string& line)
{
    auto pos = line.find('!');
    return (pos == std::string::npos) ? line : line.substr(0, pos);
}

// Parse species side: tokens like "2H + O = ..." → species → coefficient map.
std::map<std::string, double> parse_side(const std::string& side)
{
    std::map<std::string, double> out;
    std::string s = side;
    // Split on '+' but not inside parentheses (Lindemann (+M)).
    std::vector<std::string> terms;
    int depth = 0;
    std::string cur;
    for (char c : s) {
        if (c == '(')
            ++depth;
        if (c == ')')
            --depth;
        if (c == '+' && depth == 0) {
            terms.push_back(trim(cur));
            cur.clear();
        } else
            cur.push_back(c);
    }
    if (!cur.empty())
        terms.push_back(trim(cur));
    for (auto& t : terms) {
        if (t.empty() || t == "M" || t == "(+M)" || t == "+M")
            continue;
        // Leading numeric stoichiometry: "2HO2", "1.5O2", "H2"
        double coef = 1.0;
        std::size_t i = 0;
        while (i < t.size() && (std::isdigit(static_cast<unsigned char>(t[i])) || t[i] == '.'))
            ++i;
        if (i > 0 && i < t.size()) {
            coef = std::strtod(t.substr(0, i).c_str(), nullptr);
            t = t.substr(i);
        }
        out[t] += coef;
    }
    return out;
}

double energy_to_si(double Ea, EnergyUnit u)
{
    switch (u) {
    case EnergyUnit::CalPerMole:
        return Ea * 4.184;
    case EnergyUnit::KCalPerMole:
        return Ea * 4184.0;
    case EnergyUnit::JoulesPerMole:
        return Ea;
    case EnergyUnit::KJoulesPerMole:
        return Ea * 1000.0;
    case EnergyUnit::Kelvins:
        return Ea * 8.314462618;
    }
    return Ea;
}

} // namespace

ChemkinMechanism ChemkinParser::parse_file(const std::string& path) const
{
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("ChemkinParser: cannot open file: " + path);
    std::stringstream ss;
    ss << f.rdbuf();
    return parse_string(ss.str());
}

ChemkinMechanism ChemkinParser::parse_string(const std::string& text) const
{
    ChemkinMechanism M;
    std::istringstream iss(text);
    std::string raw;
    enum class Section
    {
        None,
        Elements,
        Species,
        Thermo,
        Reactions
    };
    Section sec = Section::None;
    std::string lastEquation;

    while (std::getline(iss, raw)) {
        std::string line = trim(strip_comment(raw));
        if (line.empty())
            continue;
        const std::string U = upper(line);

        if (U.find("ELEMENTS") == 0 || U.find("ELEM") == 0) {
            sec = Section::Elements;
            continue;
        }
        if (U.find("SPECIES") == 0 || U.find("SPEC") == 0) {
            sec = Section::Species;
            continue;
        }
        if (U.find("THERMO") == 0) {
            sec = Section::Thermo;
            continue;
        }
        if (U.find("REACTIONS") == 0) {
            sec = Section::Reactions;
            // Optional units header.
            auto toks = split_ws(U);
            for (std::size_t i = 1; i < toks.size(); ++i) {
                if (toks[i] == "CAL/MOLE")
                    M.energyUnit = EnergyUnit::CalPerMole;
                else if (toks[i] == "KCAL/MOLE")
                    M.energyUnit = EnergyUnit::KCalPerMole;
                else if (toks[i] == "JOULES/MOLE")
                    M.energyUnit = EnergyUnit::JoulesPerMole;
                else if (toks[i] == "KJOULES/MOLE")
                    M.energyUnit = EnergyUnit::KJoulesPerMole;
                else if (toks[i] == "KELVINS")
                    M.energyUnit = EnergyUnit::Kelvins;
                else if (toks[i] == "MOLES")
                    M.quantityUnit = QuantityUnit::Moles;
                else if (toks[i] == "MOLECULES")
                    M.quantityUnit = QuantityUnit::Molecules;
            }
            continue;
        }
        if (U == "END") {
            sec = Section::None;
            continue;
        }

        switch (sec) {
        case Section::Elements: {
            for (auto& t : split_ws(line))
                M.elements.push_back({t, atomic_weight(t)});
            break;
        }
        case Section::Species: {
            for (auto& t : split_ws(line)) {
                ChemkinSpecies sp;
                sp.name = t;
                M.species.push_back(std::move(sp));
            }
            break;
        }
        case Section::Thermo: {
            // Minimal NASA-7 detection: collect the species name on
            // 80-column record lines.  Full polynomial parsing is the
            // ThermoNasaParser's job; here we just flag which species
            // have a thermo block to enable cross-reference.
            if (line.size() >= 1 && std::isalpha(static_cast<unsigned char>(line[0]))) {
                auto tk = split_ws(line);
                if (!tk.empty()) {
                    for (auto& sp : M.species)
                        if (sp.name == tk[0])
                            sp.hasThermo = true;
                }
            }
            break;
        }
        case Section::Reactions: {
            const std::string Uline = upper(line);
            if (Uline.find("DUPLICATE") == 0 || Uline.find("DUP") == 0) {
                if (!M.reactions.empty())
                    M.reactions.back().isDuplicate = true;
                break;
            }
            if (Uline.find("LOW") == 0 && !M.reactions.empty()) {
                auto lpar = line.find('/');
                auto rpar = line.rfind('/');
                if (lpar != std::string::npos && rpar > lpar) {
                    auto v = split_ws(line.substr(lpar + 1, rpar - lpar - 1));
                    if (v.size() >= 3) {
                        M.reactions.back().isLindemann = true;
                        M.reactions.back().low.A = std::strtod(v[0].c_str(), nullptr);
                        M.reactions.back().low.beta = std::strtod(v[1].c_str(), nullptr);
                        M.reactions.back().low.Ea =
                            energy_to_si(std::strtod(v[2].c_str(), nullptr), M.energyUnit);
                    }
                }
                break;
            }
            if (Uline.find("TROE") == 0 && !M.reactions.empty()) {
                auto lpar = line.find('/');
                auto rpar = line.rfind('/');
                if (lpar != std::string::npos && rpar > lpar) {
                    auto v = split_ws(line.substr(lpar + 1, rpar - lpar - 1));
                    auto& tr = M.reactions.back().troe;
                    tr.enabled = true;
                    if (v.size() >= 1)
                        tr.a = std::strtod(v[0].c_str(), nullptr);
                    if (v.size() >= 2)
                        tr.T3 = std::strtod(v[1].c_str(), nullptr);
                    if (v.size() >= 3)
                        tr.T1 = std::strtod(v[2].c_str(), nullptr);
                    if (v.size() >= 4)
                        tr.T2 = std::strtod(v[3].c_str(), nullptr);
                }
                break;
            }
            if (Uline.find("REV") == 0 && !M.reactions.empty()) {
                auto lpar = line.find('/');
                auto rpar = line.rfind('/');
                if (lpar != std::string::npos && rpar > lpar) {
                    auto v = split_ws(line.substr(lpar + 1, rpar - lpar - 1));
                    if (v.size() >= 3) {
                        M.reactions.back().hasReverse = true;
                        M.reactions.back().rev.A = std::strtod(v[0].c_str(), nullptr);
                        M.reactions.back().rev.beta = std::strtod(v[1].c_str(), nullptr);
                        M.reactions.back().rev.Ea =
                            energy_to_si(std::strtod(v[2].c_str(), nullptr), M.energyUnit);
                    }
                }
                break;
            }
            // Third-body enhancement lines look like "H2O/16.0/ CO2/3.8/".
            if (line.find('/') != std::string::npos && Uline.find('=') == std::string::npos
                && !M.reactions.empty()) {
                std::string s = line;
                while (true) {
                    auto p1 = s.find('/');
                    if (p1 == std::string::npos)
                        break;
                    auto p2 = s.find('/', p1 + 1);
                    if (p2 == std::string::npos)
                        break;
                    std::string sp = trim(s.substr(0, p1));
                    double val = std::strtod(s.substr(p1 + 1, p2 - p1 - 1).c_str(), nullptr);
                    if (!sp.empty())
                        M.reactions.back().thirdBodyEff[sp] = val;
                    s = s.substr(p2 + 1);
                }
                break;
            }
            // Reaction line: "A + B <=> C + D    A  beta  Ea"
            if (Uline.find('=') != std::string::npos) {
                // Direction.
                ReactionDirection dir = ReactionDirection::Forward;
                std::string eq = line;
                std::string lhs, rhs;
                auto split_at = [&](const std::string& sep) {
                    auto p = eq.find(sep);
                    if (p != std::string::npos) {
                        lhs = eq.substr(0, p);
                        rhs = eq.substr(p + sep.size());
                        return true;
                    }
                    return false;
                };
                if (split_at("<=>"))
                    dir = ReactionDirection::Reversible;
                else if (split_at("=>"))
                    dir = ReactionDirection::Forward;
                else if (split_at("="))
                    dir = ReactionDirection::Reversible;
                else
                    break;
                // Trailing 3 numeric tokens are A beta Ea.
                auto toks = split_ws(rhs);
                if (toks.size() < 3)
                    break;
                ChemkinReaction r;
                r.direction = dir;
                r.fwd.Ea =
                    energy_to_si(std::strtod(toks[toks.size() - 1].c_str(), nullptr), M.energyUnit);
                r.fwd.beta = std::strtod(toks[toks.size() - 2].c_str(), nullptr);
                r.fwd.A = std::strtod(toks[toks.size() - 3].c_str(), nullptr);
                std::string rhsSpec;
                for (std::size_t i = 0; i + 3 < toks.size(); ++i) {
                    rhsSpec += toks[i];
                    rhsSpec += ' ';
                }
                r.reactants = parse_side(lhs);
                r.products = parse_side(rhsSpec);
                std::ostringstream eqOut;
                eqOut << trim(lhs) << (dir == ReactionDirection::Reversible ? " <=> " : " => ")
                      << trim(rhsSpec);
                r.equation = eqOut.str();
                lastEquation = r.equation;
                M.reactions.push_back(std::move(r));
            }
            break;
        }
        case Section::None:
            break;
        }
    }
    SIMALL_LOG_INFO("Combustion",
                    "Chemkin mechanism parsed: ",
                    M.elements.size(),
                    " elements, ",
                    M.species.size(),
                    " species, ",
                    M.reactions.size(),
                    " reactions");
    return M;
}

} // namespace simall::combustion
