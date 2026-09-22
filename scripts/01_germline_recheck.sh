#!/usr/bin/env bash
# =============================================================================
# 01_germline_recheck.sh — re-annotate the germline (blood) exome against
# CURRENT ClinVar + gnomAD, then hand off to 01_germline_recheck.R for the report.
#
# Why: Sema4's "no hereditary findings" was 2020 knowledge. ClinVar has
# reclassified tens of thousands of variants since. This is the one analysis
# with consequences for living relatives, so it goes first.
#
# Run on the Mac (needs conda; installs bcftools/htslib if missing):
#   bash 01_germline_recheck.sh
#
# Inputs : data/sema4/**/ *.germline.vcf   (Sema4's normal-sample calls, hg19)
# Refs   : refs/clinvar_GRCh37.vcf.gz (downloaded, ~100 MB), refs/refGene_hg19.txt.gz
# Output : results/01_germline/germline.annot.vcf.gz, *.tsv  (inputs to the R script)
# =============================================================================
set -euo pipefail

# Project root = the folder you run this from (must contain data/sema4), or set PROJ=...
PROJ="${PROJ:-$PWD}"
[[ -d "$PROJ/data/sema4" ]] || { echo "ERROR: $PROJ/data/sema4 not found — cd into the project folder first"; exit 1; }
DATA="$PROJ/data/sema4"
REF="$PROJ/refs"
OUT="$PROJ/results/01_germline"
mkdir -p "$REF" "$OUT"

CLINVAR_URL="${CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz}"
REFGENE_URL="${REFGENE_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz}"
# gnomAD v2.1.1 exomes (GRCh37). Read REMOTELY by region via htslib — no 58 GB download.
GNOMAD_URL="${GNOMAD_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz}"

log(){ echo "[$(date '+%H:%M:%S')] $*"; }

# ---------- tools ----------
# bcftools/htslib: Homebrew if available, otherwise a dedicated conda env (never touch base).
ENVDIR="$(conda info --base 2>/dev/null || echo /nonexistent)/envs/esos-bio"
[[ -d "$ENVDIR/bin" ]] && export PATH="$ENVDIR/bin:$PATH"
if ! command -v bcftools >/dev/null || ! command -v tabix >/dev/null; then
  if command -v brew >/dev/null; then
    log "installing bcftools/htslib via Homebrew"; brew install bcftools htslib
  else
    log "creating conda env 'esos-bio' with bcftools/htslib"
    conda create -y -n esos-bio -c conda-forge -c bioconda bcftools htslib
    export PATH="$ENVDIR/bin:$PATH"
  fi
fi
log "bcftools $(bcftools --version | head -1 | awk '{print $2}')"

# ---------- locate the germline VCF ----------
GVCF="${GVCF:-$(find -L "$DATA" -name '*.germline.vcf' ! -name '*mask_mark*' | head -1)}"
[[ -f "$GVCF" ]] || { log "ERROR: no *.germline.vcf under $DATA"; exit 1; }
log "germline VCF: $GVCF"
# Sema4 puts free-text "## ..." comment lines in the header; bcftools rejects them. Keep only ##key=value.
awk '!/^##/ || /^##[^ =]+=/' "$GVCF" > "$OUT/germline.input.vcf"
GVCF="$OUT/germline.input.vcf"
log "samples in VCF: $(bcftools query -l "$GVCF" | tr '\n' ' ')"
log "variant count : $(bcftools view -H "$GVCF" | wc -l)"

# Does Sema4 already carry a consequence annotation? (SnpEff ANN / VEP CSQ)
HAS_ANN=$(bcftools view -h "$GVCF" | grep -cE '^##INFO=<ID=(ANN|CSQ),' || true)
log "consequence annotation in header (ANN/CSQ): $HAS_ANN"

# ---------- chromosome naming: ClinVar/gnomAD use '1', not 'chr1' ----------
FIRST_CHR=$(awk '!/^#/{print $1; exit}' "$GVCF")   # no pipe: `| head -1` SIGPIPEs bcftools under pipefail
if [[ "$FIRST_CHR" == chr* ]]; then
  log "VCF uses 'chr' prefix -> stripping to match ClinVar/gnomAD"
  for c in $(seq 1 22) X Y M; do echo "chr$c $c"; done > "$OUT/chr_strip.txt"
  echo "chrM MT" >> "$OUT/chr_strip.txt"
  bcftools annotate --rename-chrs "$OUT/chr_strip.txt" "$GVCF" -Oz -o "$OUT/germline.raw.vcf.gz"
else
  bcftools view "$GVCF" -Oz -o "$OUT/germline.raw.vcf.gz"
fi
tabix -f "$OUT/germline.raw.vcf.gz"

# ---------- split multi-allelics so REF/ALT match ClinVar exactly ----------
bcftools norm -m -both "$OUT/germline.raw.vcf.gz" -Oz -o "$OUT/germline.norm.vcf.gz" 2>/dev/null
tabix -f "$OUT/germline.norm.vcf.gz"

# ---------- references ----------
if [[ ! -s "$REF/clinvar_GRCh37.vcf.gz" ]]; then
  log "downloading ClinVar (GRCh37)"
  curl -sSL -o "$REF/clinvar_GRCh37.vcf.gz" "$CLINVAR_URL"
  curl -sSL -o "$REF/clinvar_GRCh37.vcf.gz.tbi" "$CLINVAR_URL.tbi" || tabix -f "$REF/clinvar_GRCh37.vcf.gz"
fi
[[ -s "$REF/clinvar_GRCh37.vcf.gz.tbi" ]] || tabix -f "$REF/clinvar_GRCh37.vcf.gz"
log "ClinVar file date: $(bcftools view -h "$REF/clinvar_GRCh37.vcf.gz" | grep -m1 fileDate || echo unknown)"

if [[ ! -s "$REF/refGene_hg19.txt.gz" ]]; then
  log "downloading UCSC refGene (hg19) for gene coordinates"
  curl -sSL -o "$REF/refGene_hg19.txt.gz" "$REFGENE_URL"
fi

# ---------- cancer-predisposition panel (ACMG SF v3.2 cancer genes + sarcoma/OS-relevant) ----------
cat > "$OUT/panel_genes.txt" <<'EOF'
TP53
RB1
RECQL4
BLM
WRN
CDKN2A
NF1
BRCA1
BRCA2
PALB2
ATM
CHEK2
BARD1
BRIP1
RAD51C
RAD51D
NBN
MLH1
MSH2
MSH6
PMS2
EPCAM
APC
MUTYH
POLE
POLD1
PTEN
STK11
CDH1
BMPR1A
SMAD4
MEN1
RET
VHL
NF2
TSC1
TSC2
WT1
DICER1
POT1
SDHB
SDHC
SDHD
SDHAF2
MAX
TMEM127
FANCA
FANCC
FANCD2
CDK4
MDM2
EOF

# BED of panel gene bodies (+/- 20 bp), refGene name2 = gene symbol, chrom w/o 'chr'
zcat < "$REF/refGene_hg19.txt.gz" | awk -v OFS='\t' 'NR==FNR{g[$1]=1; next}
  ($13 in g) && $3 ~ /^chr([0-9]+|X|Y)$/ {c=$3; sub(/^chr/,"",c); s=$5-20; if(s<0)s=0; print c, s, $6+20, $13}' \
  "$OUT/panel_genes.txt" - | sort -k1,1V -k2,2n | bedtools merge -c 4 -o distinct -i - 2>/dev/null > "$OUT/panel.bed" \
  || zcat < "$REF/refGene_hg19.txt.gz" | awk -v OFS='\t' 'NR==FNR{g[$1]=1; next}
  ($13 in g) && $3 ~ /^chr([0-9]+|X|Y)$/ {c=$3; sub(/^chr/,"",c); s=$5-20; if(s<0)s=0; print c, s, $6+20, $13}' \
  "$OUT/panel_genes.txt" - | sort -k1,1V -k2,2n -u > "$OUT/panel.bed"
log "panel BED: $(wc -l < "$OUT/panel.bed") intervals, $(cut -f4 "$OUT/panel.bed" | tr ',' '\n' | sort -u | wc -l) genes found"

# ---------- annotate with ClinVar (exact CHROM/POS/REF/ALT match) ----------
log "annotating with ClinVar"
bcftools annotate -a "$REF/clinvar_GRCh37.vcf.gz" \
  -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/ALLELEID \
  "$OUT/germline.norm.vcf.gz" -Oz -o "$OUT/germline.clinvar.vcf.gz"
tabix -f "$OUT/germline.clinvar.vcf.gz"

# ---------- gnomAD allele frequency, panel regions only, streamed remotely ----------
# Only the panel (~50 genes) is fetched; if the remote read fails we carry on without AF.
log "annotating panel regions with gnomAD v2.1.1 exome AF (remote; may take a few minutes)"
GNOMAD_OK=0
if bcftools view -R "$OUT/panel.bed" "$GNOMAD_URL" -Oz -o "$OUT/gnomad_panel.vcf.gz" 2>"$OUT/gnomad.err"; then
  tabix -f "$OUT/gnomad_panel.vcf.gz"
  bcftools annotate -a "$OUT/gnomad_panel.vcf.gz" \
    -c 'INFO/GNOMAD_AF:=INFO/AF,INFO/GNOMAD_AF_POPMAX:=INFO/AF_popmax,INFO/GNOMAD_NHOMALT:=INFO/nhomalt' \
    "$OUT/germline.clinvar.vcf.gz" -Oz -o "$OUT/germline.annot.vcf.gz" && GNOMAD_OK=1
else
  log "WARN: gnomAD remote read failed ($(head -1 "$OUT/gnomad.err")). Continuing without AF."
fi
(( GNOMAD_OK )) || cp "$OUT/germline.clinvar.vcf.gz" "$OUT/germline.annot.vcf.gz"
tabix -f "$OUT/germline.annot.vcf.gz"

# ---------- VEP: gene / consequence / HGVS / exome-wide AF (optional; install: run_ec2.sh 00) ----------
VEP_ENV="${VEP_ENV:-/data/envs/vep}"; VEP_CACHE="${VEP_CACHE:-$REF/vep}"
FINAL="$OUT/germline.annot.vcf.gz"; HAS_VEP=0
if [[ -x "$VEP_ENV/bin/vep" && -f "$VEP_CACHE/.install_done" ]]; then
  log "annotating consequences with VEP (--pick: one transcript per variant, canonical preferred)"
  PATH="$VEP_ENV/bin:$PATH" vep -i "$OUT/germline.annot.vcf.gz" -o "$OUT/germline.vep.vcf.gz" \
    --offline --cache --dir_cache "$VEP_CACHE" --assembly GRCh37 --vcf --compress_output bgzip \
    --pick --symbol --canonical --hgvs --numbers --sift b --polyphen b --af_gnomade --max_af \
    --fork "$(nproc)" --no_stats --force_overwrite
  tabix -f "$OUT/germline.vep.vcf.gz"
  FINAL="$OUT/germline.vep.vcf.gz"; HAS_VEP=1
else
  log "VEP not installed ($VEP_ENV) — skipping consequence annotation"
fi

# ---------- flat tables for R ----------
FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%DP]\t[%AD]\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/GENEINFO\t%INFO/GNOMAD_AF\t%INFO/GNOMAD_AF_POPMAX\t%INFO/GNOMAD_NHOMALT'
HDR='CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tAD\tCLNSIG\tCLNSIGCONF\tCLNREVSTAT\tCLNDN\tGENEINFO\tGNOMAD_AF\tGNOMAD_AF_POPMAX\tGNOMAD_NHOMALT'
if (( HAS_ANN )); then FMT="$FMT\t%INFO/ANN"; HDR="$HDR\tANN"; fi
# Sema4 writes GT as a bare "1"; zygosity is in INFO/CALL_TYPE and their call QC in INFO/RED_FLAGS
for tag in CALL_TYPE RED_FLAGS; do
  if bcftools view -h "$GVCF" | grep -q "^##INFO=<ID=$tag,"; then FMT="$FMT\t%INFO/$tag"; HDR="$HDR\t$tag"; fi
done
QUERY=(bcftools query)
if (( HAS_VEP )); then   # +split-vep exposes CSQ subfields as %TAGs in the same query syntax
  QUERY=(bcftools +split-vep)
  for tag in SYMBOL Consequence IMPACT HGVSc HGVSp EXON SIFT PolyPhen gnomADe_AF MAX_AF; do FMT="$FMT\t%$tag"; HDR="$HDR\t$tag"; done
fi

# (a) every variant, genome-wide, that ClinVar knows anything about
{ printf "$HDR\n"; "${QUERY[@]}" -i 'INFO/CLNSIG!=""' -f "$FMT\n" "$FINAL"; } > "$OUT/exome_clinvar_hits.tsv"
# (b) every variant inside the panel genes, ClinVar or not
{ printf "$HDR\n"; "${QUERY[@]}" -R "$OUT/panel.bed" -f "$FMT\n" "$FINAL"; } > "$OUT/panel_all_variants.tsv"
# (c) exome-wide rare protein-truncating / splice variants ClinVar has NOT classified (needs VEP)
if (( HAS_VEP )); then
  { printf "$HDR\n"; "${QUERY[@]}" -i 'INFO/CLNSIG="."' -f "$FMT\n" "$FINAL" | awk -F'\t' -v c="$(echo -e "$HDR" | tr '\t' '\n' | grep -n '^IMPACT$' | cut -d: -f1)" '$c=="HIGH"'; } > "$OUT/exome_lof_unclassified.tsv"
  log "HIGH-impact variants not in ClinVar     : $(( $(wc -l < "$OUT/exome_lof_unclassified.tsv") - 1 ))"
fi

log "ClinVar-annotated variants (exome-wide): $(( $(wc -l < "$OUT/exome_clinvar_hits.tsv") - 1 ))"
log "variants in panel genes               : $(( $(wc -l < "$OUT/panel_all_variants.tsv") - 1 ))"
log "DONE -> now run:  Rscript 01_germline_recheck.R"
