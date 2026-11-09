# Downloads the SketchTree-C genomes from NCBI into data\
# Run from the project root (the folder that contains src, data, results):
#   powershell -ExecutionPolicy Bypass -File .\download_genomes.ps1

$genomes = [ordered]@{
    "lambda"       = "NC_001416"
    "sars2"        = "NC_045512"
    "hiv1"         = "NC_001802"
    "hbv"          = "NC_003977"
    "ebola"        = "NC_002549"
    "zika"         = "NC_012532"
    "human_mt"     = "NC_012920"
    "chimp_mt"     = "NC_001643"
    "gorilla_mt"   = "NC_001645"
    "mouse_mt"     = "NC_005089"
    "zebrafish_mt" = "NC_002333"
    "yeast_mt"     = "NC_001224"
    "mgen"         = "NC_000908"
    "ecoli"        = "NC_000913"
    "shigella"     = "NC_004337"
    "salmonella"   = "NC_003197"
    "bsub"         = "NC_000964"
    "banthracis"   = "NC_003997"
    "saureus"      = "NC_007795"
    "mtb"          = "NC_000962"
    "paeruginosa"  = "NC_002516"
    "hpylori"      = "NC_000915"
    "scoelicolor"  = "NC_003888"
    "mjannaschii"  = "NC_000909"
}

New-Item -ItemType Directory -Force -Path data | Out-Null
$failed = @()

foreach ($name in $genomes.Keys) {
    $acc = $genomes[$name]
    $out = "data\$name.fna"

    if ((Test-Path $out) -and ((Get-Item $out).Length -gt 0)) {
        Write-Host "$name already exists, skipping"
        continue
    }

    $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$acc&rettype=fasta&retmode=text"
    $ok = $false

    for ($try = 1; ($try -le 3) -and (-not $ok); $try++) {
        try {
            Invoke-WebRequest -Uri $url -OutFile $out -UseBasicParsing -ErrorAction Stop
            $first = Get-Content $out -TotalCount 1
            if ($first -and $first.StartsWith(">")) {
                $ok = $true
                Write-Host "$name  ->  $first"
            }
            else {
                Remove-Item $out -ErrorAction SilentlyContinue
                throw "NCBI did not return a FASTA record"
            }
        }
        catch {
            Write-Host "  attempt $try failed for $name ($acc): $($_.Exception.Message)" -ForegroundColor Yellow
            Start-Sleep -Seconds 3
        }
    }

    if (-not $ok) { $failed += "$name ($acc)" }
    Start-Sleep -Milliseconds 700   # NCBI allows about 3 requests per second
}

Write-Host ""
if ($failed.Count -gt 0) {
    Write-Host "These did not download:" -ForegroundColor Red
    $failed | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    Write-Host "Run the script again to retry only these."
}
else {
    Write-Host "All genomes downloaded. Check that each header above matches the organism."
}