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
    # CI-friendly mode: per-file --dry-run, collect names of files that would
    # change, print a summary list (NOT the full diffs -- the diff output for
    # this codebase exceeds 10 MB, which floods CI logs unhelpfully).
    $needs = New-Object System.Collections.Generic.List[string]
    foreach ($f in $files) {
        & clang-format --dry-run -Werror -style=file -- $f.FullName 2>$null
        if ($LASTEXITCODE -ne 0) { $needs.Add((Resolve-Path -Relative $f.FullName)) }
    }
    if ($needs.Count -eq 0) {
        Write-Host "All $($files.Count) files conform to .clang-format." -ForegroundColor Green
        exit 0
    }
    Write-Host "$($needs.Count) of $($files.Count) files would be reformatted by clang-format:" -ForegroundColor Yellow
    $needs | ForEach-Object { Write-Host "  $_" }
    Write-Host ""
    Write-Host "Run 'pwsh scripts/format.ps1' locally to fix." -ForegroundColor Yellow
    exit 1
}
else {
    & clang-format -i -style=file @($files.FullName)
    Write-Host "Formatted $($files.Count) files." -ForegroundColor Green
}
