<#
.SYNOPSIS
    Download and unpack the Inter font (OFL) into resources/fonts/inter/.

.PARAMETER Version
    Inter release tag. Default 4.0.
#>
[CmdletBinding()]
param([string]$Version = '4.0')

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
$dest = Join-Path $repo 'resources/fonts/inter'

if (Test-Path (Join-Path $dest 'Inter-Regular.ttf')) {
    Write-Host "Inter already installed." -ForegroundColor Green
    return
}

New-Item -ItemType Directory -Force -Path $dest | Out-Null
$url = "https://github.com/rsms/inter/releases/download/v$Version/Inter-$Version.zip"
$zip = Join-Path $env:TEMP "inter-$Version.zip"

Write-Host "Fetching $url"
Invoke-WebRequest -Uri $url -OutFile $zip

Expand-Archive -Force -Path $zip -DestinationPath $dest
Remove-Item $zip
Write-Host "Inter installed to $dest" -ForegroundColor Green
