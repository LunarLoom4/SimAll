# SimAll Beta — Sample Data Directory

This directory ships small, license-clean reference geometries and grids
used by the example applications and the regression suite.

Every file has a corresponding entry in `manifest.json` recording its
origin, license, intended use, and SHA-256 checksum.

## Layout

```
data/
  cad/           STEP / IGES / BREP exact geometry
  stl/           ASCII / binary STL triangulated surfaces
  mesh/          CGNS / Fluent / Plot3D pre-built grids
  fields/        reference solutions for verification (HDF5)
  manifest.json  registry of every shipped asset
```

## License

All bundled assets are either (a) authored in-house for SimAll Beta and
released under the project license, or (b) imported from public-domain
sources (NIST, NASA TM reports, Wikimedia 3D commons). Each entry in
`manifest.json` carries the exact license string.

## How to add a new asset

1. Drop the file under the appropriate subdirectory.
2. Compute its SHA-256:  `Get-FileHash <path> -Algorithm SHA256`.
3. Append a record to `manifest.json` with `origin`, `license`,
   `purpose`, `bytes`, `sha256`.
4. Run `scripts/verify_data.ps1` to confirm the manifest matches.

Large reference grids (>5 MB) MUST NOT be committed directly. Add a
download URL to `manifest.json` and let `scripts/fetch_data.ps1` pull
them on first build.
