$ErrorActionPreference = 'Stop'
$module = Split-Path -Parent $MyInvocation.MyCommand.Path
$rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
if (-not $rscript) {
    $rscript = 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe'
}
if (-not (Test-Path -LiteralPath $rscript)) {
    throw 'Rscript is required; edit $rscript in run_pipeline.ps1 for your installation.'
}
& python (Join-Path $module 'scripts\00_select_paper_biocarta.py')
if ($LASTEXITCODE -ne 0) { throw 'Paper BioCarta selection failed.' }
& $rscript (Join-Path $module 'scripts\01_export_motrpac_biocarta.R')
if ($LASTEXITCODE -ne 0) { throw 'BioCarta export failed.' }
& python (Join-Path $module 'scripts\02_map_paper_to_motrpac.py')
if ($LASTEXITCODE -ne 0) { throw 'Paper/MoTrPAC mapping failed.' }
& python (Join-Path $module 'scripts\03_plot_biocarta_bubbles.py')
if ($LASTEXITCODE -ne 0) { throw 'BioCarta figure generation failed.' }
