#!/usr/bin/env bash
# =============================================================================
# 06_hla.sh — HLA class I typing (OptiType) from normal DNA, tumor DNA and tumor RNA (EC2; ~30-60 min).
#
#   bash scripts/run_ec2.sh 06
#
# Why: HLA type is needed for any neoantigen work and for the write-up; three independent
# samples typed separately is the built-in consistency check. HLA LOH at segment level is
# already answered by FACETS (6p21.3: 3 total / 1 minor -> no LOH); this does not re-do that.
# Immune deconvolution (quanTIseq/MCP-counter) is the Mac part (06_immune.R), from genes.results.
#
# Output: results/06_hla/hla_<sample>_result.tsv (A/B/C 4-digit + objective), hla_summary.tsv
# =============================================================================
set -euo pipefail
PROJ="${PROJ:-$PWD}"; DATA="$PROJ/data/sema4"; OUT="$PROJ/results/06_hla"; mkdir -p "$OUT"
OPTI_ENV="${OPTI_ENV:-/data/envs/optitype}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/data/tmp/conda_pkgs}"
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

if [[ ! -x "$OPTI_ENV/bin/optitype" ]]; then
  log "creating conda env $OPTI_ENV (optitype)"
  mamba create -y -q -p "$OPTI_ENV" -c conda-forge -c bioconda optitype samtools
fi
export PATH="$OPTI_ENV/bin:$PATH"   # the wrapper looks up razers3 / glpsol on PATH
OPTI="$OPTI_ENV/bin/optitype"; SAM="$OPTI_ENV/bin/samtools"   # bioconda ships an "optitype" wrapper, not OptiTypePipeline.py

NBAM=$(find -L "$DATA" -name 'ISM556054-2.dedup.recal.bam' | head -1)
TBAM=$(find -L "$DATA" -name 'ISM556046-2.dedup.recal.bam' | head -1)
RBAM=$(find -L "$DATA" -name 'ISM563041-2.star.sorted.bam' | head -1)

# Reads from the extended MHC (hg19 chr6:28.5-33.5 Mb) plus unmapped, as paired FASTQ. Chrom naming per BAM.
extract(){ # <bam> <prefix>
  local bam=$1 pre=$2 chr
  chr=$("$SAM" view -H "$bam" | awk '$1=="@SQ"{sub("SN:","",$2); if($2=="6"||$2=="chr6"){print $2; exit}}')
  "$SAM" view -@4 -u "$bam" "$chr:28500000-33500000" > "$OUT/$pre.mhc.bam"
  "$SAM" view -@4 -u -f4 "$bam" > "$OUT/$pre.unmapped.bam" || true
  "$SAM" merge -@4 -f -u "$OUT/$pre.merged.bam" "$OUT/$pre.mhc.bam" "$OUT/$pre.unmapped.bam"
  "$SAM" sort -@4 -n -o "$OUT/$pre.qn.bam" "$OUT/$pre.merged.bam"
  "$SAM" fastq -@4 -1 "$OUT/$pre.R1.fq" -2 "$OUT/$pre.R2.fq" -0 /dev/null -s /dev/null -n "$OUT/$pre.qn.bam"
  rm -f "$OUT/$pre".{mhc,unmapped,merged,qn}.bam
  log "$pre: $(( $(wc -l < "$OUT/$pre.R1.fq") / 4 )) read pairs"
}

for spec in "normal_dna:$NBAM:dna" "tumor_dna:$TBAM:dna" "tumor_rna:$RBAM:rna"; do
  IFS=: read -r name bam kind <<< "$spec"
  [[ -f "$bam" ]] || { log "WARN: $name BAM missing, skipping"; continue; }
  if [[ -s "$OUT/hla_${name}_result.tsv" ]]; then log "$name: done already"; continue; fi
  log "$name: extracting MHC reads"; extract "$bam" "$name"
  log "$name: OptiType ($kind)"
  rm -rf "$OUT/opti_$name"; mkdir -p "$OUT/opti_$name"
  "$OPTI" run -i "$OUT/$name.R1.fq" -i "$OUT/$name.R2.fq" --$kind -o "$OUT/opti_$name"   # optitype 1.5 CLI (bioconda)
  cp "$(find "$OUT/opti_$name" -name '*_result.tsv' | head -1)" "$OUT/hla_${name}_result.tsv"
  rm -f "$OUT/$name".R[12].fq
done

{ printf 'sample\tA1\tA2\tB1\tB2\tC1\tC2\treads\tobjective\n'
  for f in "$OUT"/hla_*_result.tsv; do n=$(basename "$f" _result.tsv | sed 's/^hla_//'); awk -v n="$n" 'NR==2{OFS="\t"; print n,$2,$3,$4,$5,$6,$7,$8,$9}' "$f"; done; } > "$OUT/hla_summary.tsv"
column -t "$OUT/hla_summary.tsv"
log "DONE"
