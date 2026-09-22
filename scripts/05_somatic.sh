#!/usr/bin/env bash
# =============================================================================
# 05_somatic.sh — re-annotate Sema4's somatic calls with VEP (EC2 part; ~10 min).
#
#   bash scripts/run_ec2.sh 05
#
# Input : *.somatic.vcf   (Sema4 final somatic set: tumor + normal columns, INFO/RED_FLAGS QC flags,
#                          FILTER PASS/FAIL). The mutect2.twicefiltered.vcf.gz is the upstream raw set.
# Output: results/05_somatic/
#   somatic.vep.vcf.gz        VEP CSQ (one picked transcript), COSMIC/dbSNP ids, gnomAD AF
#   somatic_annotated.tsv     flat table for 05_somatic.R (CCF vs FACETS, signatures, GENIE comparison)
# =============================================================================
set -euo pipefail
PROJ="${PROJ:-$PWD}"; DATA="$PROJ/data/sema4"; REF="$PROJ/refs"; OUT="$PROJ/results/05_somatic"
mkdir -p "$OUT"
VEP_ENV="${VEP_ENV:-/data/envs/vep}"; VEP_CACHE="${VEP_CACHE:-$REF/vep}"
log(){ echo "[$(date '+%H:%M:%S')] $*"; }
[[ -x "$VEP_ENV/bin/vep" && -f "$VEP_CACHE/.install_done" ]] || { log "ERROR: VEP not installed (run_ec2.sh 00)"; exit 1; }

SVCF=$(find -L "$DATA" -name '*.somatic.vcf' ! -name '*mask_mark*' | head -1)
[[ -f "$SVCF" ]] || { log "ERROR: no somatic.vcf under $DATA"; exit 1; }
log "somatic VCF: $SVCF"
awk '!/^##/ || /^##[^ =]+=/' "$SVCF" > "$OUT/somatic.input.vcf"          # drop Sema4 free-text header lines
log "calls: $(grep -vc '^#' "$OUT/somatic.input.vcf")  PASS: $(grep -v '^#' "$OUT/somatic.input.vcf" | awk '$7=="PASS"' | wc -l)"
for c in $(seq 1 22) X Y M; do echo "chr$c $c"; done > "$OUT/chr_strip.txt"; echo "chrM MT" >> "$OUT/chr_strip.txt"
bcftools annotate --rename-chrs "$OUT/chr_strip.txt" "$OUT/somatic.input.vcf" -Ou \
  | bcftools norm -m -both -Oz -o "$OUT/somatic.norm.vcf.gz" 2>/dev/null
tabix -f "$OUT/somatic.norm.vcf.gz"

log "VEP"
PATH="$VEP_ENV/bin:$PATH" vep -i "$OUT/somatic.norm.vcf.gz" -o "$OUT/somatic.vep.vcf.gz" \
  --offline --cache --dir_cache "$VEP_CACHE" --assembly GRCh37 --vcf --compress_output bgzip \
  --pick --symbol --canonical --hgvs --numbers --sift b --polyphen b --af_gnomade --max_af \
  --check_existing --var_synonyms --fork "$(nproc)" --no_stats --force_overwrite
tabix -f "$OUT/somatic.vep.vcf.gz"

# flat table: tumor + normal AD/DP/AF, Sema4 flags, VEP fields
# per-sample [] block expands in VCF column order: tumor (ISM556046-2) then normal (ISM556054-2)
[[ "$(bcftools query -l "$OUT/somatic.vep.vcf.gz" | paste -sd, -)" == "ISM556046-2,ISM556054-2" ]] || { echo "unexpected sample order"; exit 1; }
HDR='CHROM\tPOS\tREF\tALT\tFILTER\tRED_FLAGS\tCALL_TYPE'
FMT='%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/RED_FLAGS\t%INFO/CALL_TYPE[\t%AD\t%DP\t%AF]'
HDR="$HDR\tT_AD\tT_DP\tT_AF\tN_AD\tN_DP\tN_AF"
for tag in SYMBOL Consequence IMPACT HGVSc HGVSp EXON SIFT PolyPhen Existing_variation gnomADe_AF MAX_AF CANONICAL; do FMT="$FMT\t%$tag"; HDR="$HDR\t$tag"; done
{ printf "$HDR\n"; bcftools +split-vep -f "$FMT\n" "$OUT/somatic.vep.vcf.gz"; } > "$OUT/somatic_annotated.tsv"
log "annotated rows: $(( $(wc -l < "$OUT/somatic_annotated.tsv") - 1 ))"
log "PASS by impact: $(awk -F'\t' 'NR>1 && $5=="PASS"{n[$16]++} END{for(k in n) printf "%s=%d ", k, n[k]}' "$OUT/somatic_annotated.tsv")"
log "DONE"
