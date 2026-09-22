#!/usr/bin/env bash
# =============================================================================
# 02b_requant_rna.sh — re-quantify the tumor RNA with the GDC recipe so it is like-for-like
# with TARGET-OS / TCGA-SARC (EC2; ~2 h: index build ~1 h, alignment ~45 min).
#
#   bash scripts/run_ec2.sh 02b
#
# Why: Sema4 quantified with RSEM on RefSeq/hg19; the Xena cohorts are STAR-2pass on
# GRCh38 + GENCODE v36 with GDC's TPM formula. Cross-pipeline single-gene comparisons were
# not defensible (NOTEBOOK §4 item 10). This does not fix library prep: the tumor is
# stranded total RNA (51% intronic bases), the cohorts are poly-A.
#
# Inputs : data/sema4/**/ISM563041-2_*_R{1,2}_001.fastq.gz   (150.8 M pairs, 2x101)
# Refs   : refs/grch38/  GENCODE GRCh38 primary assembly + v36 GTF + STAR index (built once)
# Output : results/02_expression/requant/
#   ISM563041-2.gdc.genes.tsv      gene_id, gene_name, gene_type, counts (unstranded / stranded),
#                                  tpm_unstranded, tpm_stranded   (GDC-style)
#   ISM563041-2.ReadsPerGene.out.tab, Log.final.out, SJ.out.tab, Chimeric.out.junction (for step 03)
# =============================================================================
set -euo pipefail
PROJ="${PROJ:-$PWD}"; DATA="$PROJ/data/sema4"; REF="$PROJ/refs/grch38"; OUT="$PROJ/results/02_expression/requant"
mkdir -p "$REF" "$OUT"
RNA_ENV="${RNA_ENV:-/data/envs/rna}"; export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-/data/tmp/conda_pkgs}"
NT=$(nproc); log(){ echo "[$(date '+%H:%M:%S')] $*"; }
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36

if [[ ! -x "$RNA_ENV/bin/STAR" ]]; then
  log "creating conda env $RNA_ENV (star 2.7.10, samtools)"
  mamba create -y -q -p "$RNA_ENV" -c conda-forge -c bioconda 'star=2.7.10b' samtools
fi
export PATH="$RNA_ENV/bin:$PATH"
log "STAR $(STAR --version)"

R1=$(find -L "$DATA" -name 'ISM563041-2_*_R1_001.fastq.gz' | head -1); R2="${R1/_R1_/_R2_}"
[[ -f "$R1" && -f "$R2" ]] || { log "ERROR: tumor RNA FASTQs not found"; exit 1; }
log "FASTQ: $R1"

# ---------- references ----------
[[ -s "$REF/GRCh38.primary_assembly.genome.fa" ]] || { log "downloading GRCh38 primary assembly (GENCODE)";
  curl -sSL "$GENCODE/GRCh38.primary_assembly.genome.fa.gz" | gunzip -c > "$REF/GRCh38.primary_assembly.genome.fa"; }
[[ -s "$REF/gencode.v36.annotation.gtf" ]] || { log "downloading GENCODE v36 GTF";
  curl -sSL "$GENCODE/gencode.v36.annotation.gtf.gz" | gunzip -c > "$REF/gencode.v36.annotation.gtf"; }
if [[ ! -s "$REF/star_index/SA" ]]; then
  log "building STAR index (sjdbOverhang 100; ~1 h, ~32 GB RAM)"
  mkdir -p "$REF/star_index"
  STAR --runMode genomeGenerate --runThreadN "$NT" --genomeDir "$REF/star_index" \
       --genomeFastaFiles "$REF/GRCh38.primary_assembly.genome.fa" --sjdbGTFfile "$REF/gencode.v36.annotation.gtf" \
       --sjdbOverhang 100 --outTmpDir /data/tmp/star_index_tmp
fi

# ---------- STAR 2-pass, GDC parameters, gene counts + chimeric junctions ----------
if [[ ! -s "$OUT/ISM563041-2.ReadsPerGene.out.tab" ]]; then
  log "STAR alignment (2-pass) of 150.8 M pairs"
  rm -rf /data/tmp/star_run; mkdir -p /data/tmp/star_run
  STAR --runThreadN "$NT" --genomeDir "$REF/star_index" --readFilesIn "$R1" "$R2" --readFilesCommand zcat \
       --outFileNamePrefix "$OUT/ISM563041-2." --outTmpDir /data/tmp/star_run/tmp \
       --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 \
       --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 \
       --quantMode GeneCounts --outSAMtype BAM Unsorted --outSAMunmapped None --outSAMattributes NH HI AS nM NM \
       --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 --chimOutType Junctions \
       --chimMultimapNmax 10 --chimScoreJunctionNonGTAG -4 --alignSJstitchMismatchNmax 5 -1 5 5 \
       --limitBAMsortRAM 40000000000
fi
grep -E 'input reads|Uniquely mapped reads %|mapped to multiple loci %|chimeric' "$OUT/ISM563041-2.Log.final.out" | sed 's/^ *//'

# ---------- GDC-style TPM: counts / union-exon gene length, per strandedness column ----------
# ReadsPerGene columns: gene, unstranded, counts for read1 strand, counts for read2 strand (= reverse-stranded libraries)
log "computing TPM (union exon lengths from GTF)"
python3 - "$REF/gencode.v36.annotation.gtf" "$OUT/ISM563041-2.ReadsPerGene.out.tab" "$OUT/ISM563041-2.gdc.genes.tsv" <<'EOF'
import sys, re, collections
gtf, counts, outp = sys.argv[1:4]
exons = collections.defaultdict(list); name = {}; gtype = {}
for line in open(gtf):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if f[2] != 'exon': continue
    gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    exons[gid].append((int(f[3]), int(f[4])))
    if gid not in name:
        name[gid] = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        gtype[gid] = re.search(r'gene_type "([^"]+)"', f[8]).group(1)
length = {}
for gid, ivs in exons.items():
    ivs.sort(); tot = 0; cs, ce = ivs[0]
    for s, e in ivs[1:]:
        if s <= ce + 1: ce = max(ce, e)
        else: tot += ce - cs + 1; cs, ce = s, e
    length[gid] = tot + ce - cs + 1
rows = []
for line in open(counts):
    g, un, s1, s2 = line.rstrip('\n').split('\t')
    if g.startswith('N_'): continue
    rows.append((g, int(un), int(s1), int(s2)))
def tpm(col):
    rpk = {g: r[col] / (length[g] / 1000) for r in rows for g in [r[0]] if g in length}
    tot = sum(rpk.values()) / 1e6
    return {g: v / tot for g, v in rpk.items()}
tu, ts = tpm(1), tpm(3)      # unstranded (GDC convention for cohort TPM) and reverse-stranded (this library)
with open(outp, 'w') as o:
    o.write('gene_id\tgene_name\tgene_type\tlength\tunstranded\tstranded_second\ttpm_unstranded\ttpm_stranded\n')
    for g, un, s1, s2 in rows:
        if g in length:
            o.write(f'{g}\t{name[g]}\t{gtype[g]}\t{length[g]}\t{un}\t{s2}\t{tu[g]:.4f}\t{ts[g]:.4f}\n')
print('genes written:', sum(1 for r in rows if r[0] in length))
EOF
# strandedness sanity: for a reverse-stranded library column 4 >> column 3
awk 'NR>4{a+=$3; b+=$4} END{printf "assigned reads: read1-strand %d, read2-strand %d  (read2 >> read1 confirms reverse-stranded)\n", a, b}' "$OUT/ISM563041-2.ReadsPerGene.out.tab"
rm -f "$OUT/ISM563041-2.Aligned.out.bam"     # 15 GB; re-creatable, not needed downstream
log "DONE"
