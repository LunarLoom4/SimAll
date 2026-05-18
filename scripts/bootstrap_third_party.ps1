<#
.SYNOPSIS
    Bootstrap header-only / fetchable third-party dependencies via
    `git submodule` so the project builds offline once initialised.

.NOTES
    Networked. Re-run when third_party/submodules.cmake changes.
#>
[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo

try {
    if (-not (Test-Path .git)) {
        throw "Run from a git working tree (no .git found)."
    }
    Write-Host "Updating submodules…" -ForegroundColor Cyan
    git submodule sync --recursive
    git submodule update --init --recursive --depth 1
    Write-Host "Submodules synced." -ForegroundColor Green
}
finally {
    Pop-Location
}
