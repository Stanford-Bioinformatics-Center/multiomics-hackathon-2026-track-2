param([string]$RscriptPath = '')
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if (-not $RscriptPath) {
    $onPath = Get-Command Rscript -ErrorAction SilentlyContinue
    if ($onPath) {
        $RscriptPath = $onPath.Source
    } else {
        $installRoot = 'C:\Program Files\R'
        if (Test-Path -LiteralPath $installRoot) {
            $RscriptPath = Get-ChildItem -LiteralPath $installRoot -Directory -Filter 'R-*' |
                Sort-Object Name -Descending |
                ForEach-Object { Join-Path $_.FullName 'bin\Rscript.exe' } |
                Where-Object { Test-Path -LiteralPath $_ } |
                Select-Object -First 1
        }
    }
}
if (-not $RscriptPath -or -not (Test-Path -LiteralPath $RscriptPath)) {
    throw 'Rscript not found. Install R or pass -RscriptPath with its executable path.'
}
Push-Location -LiteralPath $root
try {
    & (Join-Path $PSScriptRoot '00_verify_inputs.ps1')
    foreach ($name in @('01_map_paper_to_rat.R', '02_prepare_rat_samples.R',
                        '03_estimate_training_course.R', '04_compare_pooling_choices.R',
                        '05_plot_training_course.R')) {
        Write-Output "Running $name"
        & $RscriptPath (Join-Path $PSScriptRoot $name)
        if ($LASTEXITCODE -ne 0) { throw "$name failed with exit code $LASTEXITCODE" }
    }
    Write-Output 'Pipeline complete. Results are in output/.'
}
finally { Pop-Location }
