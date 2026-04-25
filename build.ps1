# MatchY release build script
# Usage: .\build.ps1
# Produces: NSIS installer, MSI installer, and portable zip in C:\cargo-target\matchy\release\bundle\

$ErrorActionPreference = "Stop"

$AppDir   = "$PSScriptRoot\matchy\crates\matchy-app"
$TauriBin = "$AppDir\frontend\node_modules\.bin\tauri.cmd"
$TargetDir = "C:\cargo-target\matchy"
$BundleDir = "$TargetDir\release\bundle"

Write-Host "==> Building MatchY release..." -ForegroundColor Cyan

Set-Location $AppDir
$env:CARGO_TARGET_DIR = $TargetDir

& $TauriBin build
if ($LASTEXITCODE -ne 0) { Write-Error "tauri build failed"; exit 1 }

# Tauri has no 'portable' bundle target — create the zip manually from the compiled exe.
$ExePath   = "$TargetDir\release\matchy-app.exe"
$PortableDir = "$BundleDir\portable"
$PortableZip = "$PortableDir\MatchY_0.1.0_x64_portable.zip"

New-Item -ItemType Directory -Force -Path $PortableDir | Out-Null
Compress-Archive -Path $ExePath -DestinationPath $PortableZip -Force

Write-Host ""
Write-Host "==> Build complete. Artifacts:" -ForegroundColor Green
Get-ChildItem -Recurse -File $BundleDir | ForEach-Object {
    Write-Host "    $($_.FullName)"
}
