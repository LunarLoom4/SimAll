// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/Janaf.cpp
// =============================================================================
#include "materials/Janaf.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <istream>
#include <sstream>
#include <string>

namespace simall::materials
{

namespace
{

std::string toUpper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
    return s;
}

std::string rstrip(std::string s)
{
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
        s.pop_back();
    return s;
}

/// CHEMKIN NASA-7 uses Fortran D/E format. Parse a fixed-width 15-char block.
bool parseFortranDouble(const std::string& s, double& out)
{
    std::string t = s;
    // Fortran sometimes uses 'D' for exponent.
    for (auto& c : t)
        if (c == 'D' || c == 'd')
            c = 'E';
    try {
        out = std::stod(t);
        return true;
    } catch (...) {
        return false;
    }
}

} // namespace

bool JanafThermoParser::parseRecord(std::istream& in, JanafSpecies& out)
{
    std::string line1, line2, line3, line4;

    // Skip blank lines.
    do {
        if (!std::getline(in, line1))
            return false;
    } while (line1.find_first_not_of(" \t\r\n") == std::string::npos);

    std::string trimmed = line1;
    while (!trimmed.empty() && std::isspace(static_cast<unsigned char>(trimmed.front())))
        trimmed.erase(trimmed.begin());

    if (trimmed.substr(0, 3) == "END" || trimmed.substr(0, 3) == "end")
        return false;

    if (!std::getline(in, line2) || !std::getline(in, line3) || !std::getline(in, line4)) {
        errors_.push_back("JANAF: truncated record for " + line1);
        return false;
    }

    // Pad short lines to 80 chars to avoid substr exceptions.
    auto pad = [](std::string& s) {
        if (s.size() < 80)
            s.append(80 - s.size(), ' ');
    };
    pad(line1);
    pad(line2);
    pad(line3);
    pad(line4);

    // Line 1: cols 1-18 species name, 45-73 element composition (skipped),
    // col 79 phase, cols 45-54 Tlow, 55-64 Thigh, 65-72 Tmid.
    out.name = toUpper(rstrip(line1.substr(0, 18)));

    double tlow = 0, thigh = 0, tmid = 0;
    parseFortranDouble(line1.substr(45, 10), tlow);
    parseFortranDouble(line1.substr(55, 10), thigh);
    parseFortranDouble(line1.substr(65, 8), tmid);
    if (tlow > 0)
        out.Tmin = tlow;
    if (thigh > 0)
        out.Tmax = thigh;
    if (tmid > 0)
        out.Tmid = tmid;

    auto readCoeffs =
        [&](const std::string& src, std::size_t startCol, std::size_t count, double* dst) {
            for (std::size_t i = 0; i < count; ++i) {
                double v = 0.0;
                parseFortranDouble(src.substr(startCol + i * 15, 15), v);
                dst[i] = v;
            }
        };

    // Line 2: high-range a1..a5  (cols 1-75)
    readCoeffs(line2, 0, 5, out.highCoeffs.data());
    // Line 3: high a6,a7 then low a1,a2,a3
    readCoeffs(line3, 0, 2, out.highCoeffs.data() + 5);
    readCoeffs(line3, 30, 3, out.lowCoeffs.data());
    // Line 4: low a4..a7
    readCoeffs(line4, 0, 4, out.lowCoeffs.data() + 3);

    if (out.molecularWeight == 0.0) {
        // CHEMKIN .dat doesn't carry MW; the catalog supplies it later.
        out.molecularWeight = 0.0;
    }
    return true;
}

std::unordered_map<std::string, JanafSpecies> JanafThermoParser::parse(std::istream& in)
{
    std::unordered_map<std::string, JanafSpecies> out;

    std::string header;
    // Skip the optional "THERMO" / "THERMO ALL" line and the global
    // T-range line that follows it (if present). The convention is brittle;
    // we look one line ahead.
    auto pos = in.tellg();
    if (std::getline(in, header)) {
        const std::string upper = toUpper(header);
        if (upper.find("THERMO") == std::string::npos) {
            in.seekg(pos); // not a header — push back
        } else {
            // Eat the global temperature line (200.0 1000.0 5000.0).
            std::string tline;
            std::getline(in, tline);
        }
    }

    JanafSpecies sp;
    while (parseRecord(in, sp)) {
        if (!sp.name.empty())
            out[sp.name] = sp;
        sp = JanafSpecies{};
    }
    return out;
}

} // namespace simall::materials
