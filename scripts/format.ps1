<#
.SYNOPSIS
    Run clang-format -i over the entire source tree.

.PARAMETER Check
    If supplied, exits non-zero when files would change (CI mode).
#>
[CmdletBinding()]
param([switch]$Check)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

$tool = (Get-Command clang-format -ErrorAction SilentlyContinue)
if (-not $tool) { throw "clang-format not found in PATH." }

$paths = @('src', 'tests', 'applications', 'plugins')
$exts  = '*.cpp', '*.hpp', '*.h', '*.cc', '*.cu', '*.cuh'

$files = foreach ($p in $paths) {
    $abs = Join-Path $repo $p
    if (Test-Path $abs) {
        Get-ChildItem -Path $abs -Recurse -Include $exts -File
    }
}

if ($Check) {
    # clang-format >= 10 supports --dry-run -Werror; non-zero exit if any
    # file would be modified. This is robust to line-ending differences
    # whereas an in-process diff is not.
    & clang-format --dry-run -Werror -style=file @($files.FullName)
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    Write-Host "All formatted." -ForegroundColor Green
}
else {
    & clang-format -i -style=file @($files.FullName)
    Write-Host "Formatted $($files.Count) files." -ForegroundColor Green
}
