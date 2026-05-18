// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/ThermoNasaParser.cpp
// =============================================================================
#include "combustion/ThermoNasaParser.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace simall::combustion {

namespace {

std::string trim(std::string s) {
    auto issp = [](unsigned char c){ return std::isspace(c) != 0; };
    while (!s.empty() && issp(s.front())) s.erase(s.begin());
    while (!s.empty() && issp(s.back()))  s.pop_back();
    return s;
}

double atomic_weight(const std::string& sym) {
    static const std::map<std::string, double> tbl = {
        {"H", 1.00794},  {"D", 2.0141},   {"C", 12.0107},
        {"N", 14.0067},  {"O", 15.9994},  {"AR", 39.948},
        {"HE", 4.0026},  {"NE", 20.1797}, {"S", 32.065},
        {"F", 18.9984},  {"CL", 35.453},  {"BR", 79.904},
        {"I", 126.9045}, {"P", 30.9738},  {"SI", 28.0855}
    };
    std::string u = sym;
    std::transform(u.begin(), u.end(), u.begin(),
                   [](unsigned char c){ return std::toupper(c); });
    auto it = tbl.find(u);
    return (it == tbl.end()) ? 0.0 : it->second;
}

double sub_atof(const std::string& s, std::size_t a, std::size_t b) {
    if (a >= s.size()) return 0.0;
    b = std::min(b, s.size());
    return std::strtod(trim(s.substr(a, b-a)).c_str(), nullptr);
}

// Pick the temperature range whose [Tlow, Thigh] contains T (or clamp).
const NasaTempRange& pick_range(const NasaSpeciesThermo& sp, double T) {
    for (const auto& r : sp.ranges) {
        if (T >= r.Tlow && T <= r.Thigh) return r;
    }
    // Clamp: pick boundary range.
    if (T < sp.ranges.front().Tlow)  return sp.ranges.front();
    return sp.ranges.back();
}

}  // namespace

double ThermoNasaParser::cp(const NasaSpeciesThermo& sp, double T) {
    if (sp.ranges.empty() || sp.molarMass <= 0.0) return 0.0;
    const auto& r = pick_range(sp, T);
    double cp_R;
    if (sp.isNasa9) {
        cp_R = r.a[0]/(T*T) + r.a[1]/T + r.a[2] + r.a[3]*T
             + r.a[4]*T*T + r.a[5]*T*T*T + r.a[6]*T*T*T*T;
    } else {
        cp_R = r.a[0] + r.a[1]*T + r.a[2]*T*T + r.a[3]*T*T*T + r.a[4]*T*T*T*T;
    }
    return cp_R * Ru / (sp.molarMass * 1.0e-3);   // J/(kg·K), molarMass in g/mol
}

double ThermoNasaParser::h(const NasaSpeciesThermo& sp, double T) {
    if (sp.ranges.empty() || sp.molarMass <= 0.0) return 0.0;
    const auto& r = pick_range(sp, T);
    double h_RT;
    if (sp.isNasa9) {
        h_RT = -r.a[0]/(T*T) + r.a[1]*std::log(T)/T + r.a[2]
             + r.a[3]*T/2.0 + r.a[4]*T*T/3.0 + r.a[5]*T*T*T/4.0
             + r.a[6]*T*T*T*T/5.0 + r.a[7]/T;
    } else {
        h_RT = r.a[0] + r.a[1]*T/2.0 + r.a[2]*T*T/3.0
             + r.a[3]*T*T*T/4.0 + r.a[4]*T*T*T*T/5.0 + r.a[5]/T;
    }
    return h_RT * Ru * T / (sp.molarMass * 1.0e-3);
}

double ThermoNasaParser::s(const NasaSpeciesThermo& sp, double T) {
    if (sp.ranges.empty() || sp.molarMass <= 0.0) return 0.0;
    const auto& r = pick_range(sp, T);
    double s_R;
    if (sp.isNasa9) {
        s_R = -r.a[0]/(2.0*T*T) - r.a[1]/T + r.a[2]*std::log(T)
            + r.a[3]*T + r.a[4]*T*T/2.0 + r.a[5]*T*T*T/3.0
            + r.a[6]*T*T*T*T/4.0 + r.a[8];
    } else {
        s_R = r.a[0]*std::log(T) + r.a[1]*T + r.a[2]*T*T/2.0
            + r.a[3]*T*T*T/3.0 + r.a[4]*T*T*T*T/4.0 + r.a[6];
    }
    return s_R * Ru / (sp.molarMass * 1.0e-3);
}

std::vector<NasaSpeciesThermo>
ThermoNasaParser::parse_file(const std::string& path) const {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("ThermoNasaParser: cannot open " + path);
    std::stringstream ss; ss << f.rdbuf();
    return parse_string(ss.str());
}

std::vector<NasaSpeciesThermo>
ThermoNasaParser::parse_string(const std::string& text) const {
    std::vector<NasaSpeciesThermo> out;
    std::istringstream iss(text);
    std::vector<std::string> lines;
    std::string raw;
    while (std::getline(iss, raw)) lines.push_back(raw);

    // NASA-7 records: 4 consecutive lines, distinguished by columns 80 = '1','2','3','4'.
    auto parse_nasa7_block = [&](std::size_t i) -> std::size_t {
        if (i + 3 >= lines.size()) return i + 1;
        const auto& L1 = lines[i];
        if (L1.size() < 80 || L1[79] != '1') return i + 1;
        NasaSpeciesThermo sp; sp.isNasa9 = false;
        // Line 1 columns:
        //   1-18  : species name
        //   19-24 : date (ignored)
        //   25-44 : element composition  4 × (2 char element, 3 digit count)
        //   45    : phase (G/L/S)
        //   46-55 : Tlow
        //   56-65 : Thigh
        //   66-73 : Tcommon
        //   80    : '1'
        sp.name = trim(L1.substr(0, 18));
        for (int k = 0; k < 4; ++k) {
            const std::size_t off = 24 + k * 5;
            if (L1.size() < off + 5) break;
            std::string sym = trim(L1.substr(off, 2));
            std::string cnt = trim(L1.substr(off + 2, 3));
            if (sym.empty() || cnt.empty() || cnt == "0") continue;
            const int n = std::atoi(cnt.c_str());
            if (n != 0) sp.composition[sym] += n;
        }
        const double Tlow  = sub_atof(L1, 45, 55);
        const double Thigh = sub_atof(L1, 55, 65);
        const double Tmid  = sub_atof(L1, 65, 73);
        // Lines 2-4 give 14 coefficients in 5 × 15-column fields.
        auto read_coefs = [&](std::size_t row, double* dst, int count) {
            const std::string& Ln = lines[row];
            for (int k = 0; k < count; ++k) {
                dst[k] = sub_atof(Ln, k * 15, (k + 1) * 15);
            }
        };
        double c[15]{};
        read_coefs(i + 1, c + 0, 5);
        read_coefs(i + 2, c + 5, 5);
        read_coefs(i + 3, c + 10, 5);
        NasaTempRange hi; hi.Tlow = Tmid;  hi.Thigh = Thigh;
        NasaTempRange lo; lo.Tlow = Tlow;  lo.Thigh = Tmid;
        for (int k = 0; k < 7; ++k) hi.a[k] = c[k];
        for (int k = 0; k < 7; ++k) lo.a[k] = c[7 + k];
        sp.ranges.push_back(lo);
        sp.ranges.push_back(hi);
        for (auto& [sym, n] : sp.composition)
            sp.molarMass += n * atomic_weight(sym);
        out.push_back(std::move(sp));
        return i + 4;
    };

    for (std::size_t i = 0; i < lines.size(); ) {
        const std::string U = lines[i];
        const std::string trimmedU = trim(U);
        if (trimmedU.empty()
            || trimmedU.rfind("THERMO", 0) == 0
            || trimmedU.rfind("END",    0) == 0
            || trimmedU[0] == '!') { ++i; continue; }
        // NASA-7 detection: 80-char record with '1' at column 80.
        if (U.size() >= 80 && U[79] == '1') {
            i = parse_nasa7_block(i);
        } else {
            ++i;
        }
    }
    SIMALL_LOG_INFO("Combustion",
        "NASA thermo parsed: ", out.size(), " species records");
    return out;
}

}  // namespace simall::combustion
