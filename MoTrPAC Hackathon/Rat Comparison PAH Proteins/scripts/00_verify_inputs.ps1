$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$paper = Join-Path $root 'data\paper\pah_lower_proteins_malenfant2015.csv'
$paperPdf = Join-Path $root 'data\paper\Malenfant2015.pdf'
$ratDir = Join-Path $root 'data\rat'
$manifest = Get-Content -LiteralPath (Join-Path $ratDir 'source_manifest.json') -Raw | ConvertFrom-Json

if (-not (Test-Path -LiteralPath $paperPdf)) { throw 'Malenfant2015.pdf is missing.' }
$rows = @(Import-Csv -LiteralPath $paper)
if ($rows.Count -ne 9) { throw "Expected nine Table 2 downregulated proteins; found $($rows.Count)." }
if (@($rows.paper_symbol | Sort-Object -Unique).Count -ne 9) { throw 'Paper symbols are not unique.' }
if (@($rows.uniprot_accession | Sort-Object -Unique).Count -ne 9) { throw 'UniProt accessions are not unique.' }
foreach ($row in $rows) {
    $ratio = [double]::Parse($row.pah_to_control_ratio, [cultureinfo]::InvariantCulture)
    if ($ratio -ge 1 -or $ratio -le 0) { throw "Unexpected PAH/control ratio for $($row.paper_symbol)." }
}

foreach ($entry in $manifest.files.PSObject.Properties) {
    $path = Join-Path $ratDir $entry.Name
    if (-not (Test-Path -LiteralPath $path)) { throw "Missing rat data: $($entry.Name)." }
    $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $entry.Value) { throw "SHA-256 mismatch for $($entry.Name)." }
    Write-Output "Verified $($entry.Name)"
}
Write-Output 'Verified the nine paper rows and all rat data files.'
