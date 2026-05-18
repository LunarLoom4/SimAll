<#
.SYNOPSIS
    SimAll Beta — build driver (Windows / PowerShell 5.1+).

.PARAMETER Preset
    A CMake configure preset. Default: "default" (MSVC Release+AVX2+Tests).
    Available: default, debug, linux-release (WSL), hpc.

.PARAMETER Target
    Optional build target (e.g. "simall_solver", "simall_beta", "install").

.PARAMETER Jobs
    Parallel jobs. Default = number of logical cores.

.PARAMETER Clean
    Delete build/ before configuring.

.EXAMPLE
    ./scripts/build.ps1
.EXAMPLE
    ./scripts/build.ps1 -Preset debug -Target simall_solver
.EXAMPLE
    ./scripts/build.ps1 -Clean -Preset default
#>
[CmdletBinding()]
param(
    [string]$Preset = 'default',
    [string]$Target = '',
    [int]   $Jobs   = [Environment]::ProcessorCount,
    [switch]$Clean
)

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot

Push-Location $repo
try {
    if ($Clean -and (Test-Path build)) {
        Write-Host "Removing build/" -ForegroundColor Yellow
        Remove-Item -Recurse -Force build
    }

    Write-Host "Configuring preset '$Preset'" -ForegroundColor Cyan
    cmake --preset $Preset
    if ($LASTEXITCODE -ne 0) { throw "Configure failed." }

    $buildArgs = @('--build', "build")
    switch ($Preset) {
        'default' { $buildArgs += @('--config', 'Release') }
        'debug'   { $buildArgs += @('--config', 'Debug') }
    }
    if ($Target) { $buildArgs += @('--target', $Target) }
    $buildArgs += @('-j', "$Jobs")

    Write-Host "Building: cmake $($buildArgs -join ' ')" -ForegroundColor Cyan
    cmake @buildArgs
    if ($LASTEXITCODE -ne 0) { throw "Build failed." }

    Write-Host "OK" -ForegroundColor Green
}
finally {
    Pop-Location
}
