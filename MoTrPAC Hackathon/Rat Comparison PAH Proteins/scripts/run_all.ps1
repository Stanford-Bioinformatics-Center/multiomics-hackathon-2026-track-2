param([string]$RscriptPath = 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe')
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
if (-not (Test-Path -LiteralPath $RscriptPath)) {
    $RscriptPath = (Get-Command Rscript -ErrorAction Stop).Source
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
