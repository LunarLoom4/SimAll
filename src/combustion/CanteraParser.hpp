// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/CanteraParser.hpp
// Phase  : 11 — Cantera-compatible mechanism parser.
//
// Cantera (https://cantera.org) is the de-facto open-source chemistry
// toolkit and ships mechanisms in YAML (preferred) and the legacy CTI text
// format.  Section 11 of the master blueprint mandates Cantera-compatible
// parsing alongside the CHEMKIN parser.  This reader transcribes the YAML
// (or CTI) representation of species, NASA-7 thermodynamics, transport
// data, and Arrhenius reaction rates into the same neutral
// `ChemkinMechanism` record that the rest of the combustion pipeline
// already consumes, so downstream solvers see a single mechanism schema
// regardless of the source format.
//
// The implementation intentionally keeps the parser self-contained (no
// libyaml dependency) — only a minimalist YAML subset (block mappings,
// block sequences, scalars, comments) is supported, which covers every
// reaction mechanism shipped by Cantera as of v3.0.
// =============================================================================
#pragma once

#include "combustion/ChemkinParser.hpp"

#include <string>

namespace simall::combustion
{

enum class CanteraFormat
{
    Auto,
    Yaml,
    Cti
};

struct CanteraReadOptions
{
    CanteraFormat format = CanteraFormat::Auto; ///< sniff by extension if Auto
    bool stripUnused = false;                   ///< drop species not referenced by any reaction
    bool validateThermo = true;                 ///< sanity-check NASA-7 polynomials
};

class CanteraParser
{
public:
    /// Parse a Cantera mechanism file.  Returns the resulting
    /// `ChemkinMechanism` ready to be handed to LaminarFiniteRate /
    /// FgmFlamelet / TransportedPdf, etc.  Throws on IO or parse failure.
    ChemkinMechanism parse(const std::string& path, const CanteraReadOptions& opts = {});

    /// Parse from an in-memory string (useful for unit tests and embedded
    /// mechanism literals).
    ChemkinMechanism parseString(const std::string& text,
                                 CanteraFormat fmt = CanteraFormat::Yaml,
                                 const CanteraReadOptions& opts = {});
};

} // namespace simall::combustion
