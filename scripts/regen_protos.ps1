<#
.SYNOPSIS
    Regenerate protobuf bindings (placeholder — protobuf is not
    yet a hard dependency; the script becomes active once
    src/io/proto/ exists).
#>
[CmdletBinding()] param()

$repo = Split-Path -Parent $PSScriptRoot
$proto = Join-Path $repo 'src/io/proto'
$out   = Join-Path $repo 'src/io/proto/generated'

if (-not (Test-Path $proto)) {
    Write-Host "No proto sources yet — nothing to do." -ForegroundColor Yellow
    return
}

if (-not (Get-Command protoc -ErrorAction SilentlyContinue)) {
    throw "protoc not found in PATH."
}

New-Item -ItemType Directory -Force -Path $out | Out-Null

Get-ChildItem -Path $proto -Filter *.proto | ForEach-Object {
    Write-Host "protoc $($_.Name)"
    & protoc --cpp_out=$out --proto_path=$proto $_.FullName
}
Write-Host "Done." -ForegroundColor Green
