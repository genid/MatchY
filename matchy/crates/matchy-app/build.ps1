# MatchY build script — builds frontend + Rust backend via Tauri CLI
# Usage: .\build.ps1
# Output: C:\cargo-target\matchy\release\matchy-app.exe

$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

Write-Host "==> Building frontend..." -ForegroundColor Cyan
Set-Location "$scriptDir\frontend"
npm run build
if ($LASTEXITCODE -ne 0) { Write-Error "Frontend build failed"; exit 1 }

Write-Host "==> Building Rust backend (release, no installer)..." -ForegroundColor Cyan
Set-Location $scriptDir
& "$scriptDir\frontend\node_modules\.bin\tauri.cmd" build --no-bundle
if ($LASTEXITCODE -ne 0) { Write-Error "Rust build failed"; exit 1 }

Write-Host "==> Done! Binary: C:\cargo-target\matchy\release\matchy-app.exe" -ForegroundColor Green
