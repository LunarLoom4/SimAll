# =============================================================================
# SimAll Beta - scripts/package.ps1
# Windows-side wrapper: configure + build + cpack producing NSIS + ZIP under
# build\_packages\. Equivalent to scripts/package.sh on Linux/macOS.
# =============================================================================
[CmdletBinding()]
param(
    [string]$Preset    = "default",
    [string]$Generator = "NSIS;ZIP"
)
$ErrorActionPreference = 'Stop'

$Root      = Resolve-Path (Join-Path $PSScriptRoot "..")
$BuildDir  = Join-Path $Root  "build"
$OutDir    = Join-Path $BuildDir "_packages"

Write-Host "[package.ps1] configure preset=$Preset"
cmake --preset $Preset

Write-Host "[package.ps1] build"
cmake --build $BuildDir --config Release --parallel

Write-Host "[package.ps1] cpack ($Generator)"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
Push-Location $BuildDir
try {
    cpack -B $OutDir -G $Generator --config CPackConfig.cmake
} finally { Pop-Location }

Write-Host "[package.ps1] artefacts:"
Get-ChildItem $OutDir | Select-Object Name,Length,LastWriteTime
