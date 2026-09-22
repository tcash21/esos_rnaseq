#!/usr/bin/env bash
# =============================================================================
# 00_setup_vep.sh — one-time install of Ensembl VEP + GRCh37 offline cache on EC2.
#
#   bash scripts/run_ec2.sh 00        # ~30-60 min, dominated by the ~17 GB cache download
#
# VEP turns CHROM/POS/REF/ALT into gene, transcript, consequence (missense, stop_gained,
# splice...), HGVS c./p. names, SIFT/PolyPhen, and population AF. Steps 01 (germline) and
# 05 (somatic) use it when present and skip it when not.
#
# Own conda prefix (/data/envs/vep) so its perl/htslib never shadow the 'bio' env tools.
# =============================================================================
set -euo pipefail

PROJ="${PROJ:-$PWD}"
VEP_ENV="${VEP_ENV:-/data/envs/vep}"
VEP_CACHE="${VEP_CACHE:-$PROJ/refs/vep}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/data/tmp/conda_pkgs}"
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

if [[ ! -x "$VEP_ENV/bin/vep" ]]; then
  log "creating conda env $VEP_ENV (ensembl-vep)"
  mamba create -y -q -p "$VEP_ENV" -c conda-forge -c bioconda ensembl-vep
fi
export PATH="$VEP_ENV/bin:$PATH"
log "VEP $(vep --help 2>/dev/null | grep -m1 -oE 'ensembl-vep +: +[0-9.]+' || echo '?')"

if [[ ! -f "$VEP_CACHE/.install_done" ]]; then   # marker, not the dir: an interrupted run leaves a partial cache
  log "downloading GRCh37 cache + FASTA into $VEP_CACHE"
  rm -rf "$VEP_CACHE"; mkdir -p "$VEP_CACHE"
  vep_install -a cf -s homo_sapiens -y GRCh37 -c "$VEP_CACHE" --NO_UPDATE
  touch "$VEP_CACHE/.install_done"
fi
log "cache: $(ls -d "$VEP_CACHE"/homo_sapiens/* | tr '\n' ' ')  ($(du -sh "$VEP_CACHE" | cut -f1))"

# smoke test: TP53 R175H (hg19 17:7578406 C>T) must come back as missense in TP53
printf '17\t7578406\t.\tC\tT\t.\t.\t.\n' > /data/tmp/vep_smoke.vcf
vep -i /data/tmp/vep_smoke.vcf -o STDOUT --offline --cache --dir_cache "$VEP_CACHE" --assembly GRCh37 \
    --pick --symbol --hgvs --tab --no_stats --fields "SYMBOL,Consequence,HGVSp" 2>/dev/null | grep -v '^#'
log "DONE"
