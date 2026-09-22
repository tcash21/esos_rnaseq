#!/usr/bin/env bash
# =============================================================================
# 04_copynumber.sh — allele-specific copy number from the tumor/normal BAMs (EC2).
#
#   bash scripts/run_ec2.sh 04       # ~1-2 h, dominated by snp-pileup over the 34 GB tumor BAM
#
# Why: Sema4's GATK calls assume a diploid baseline. The segment file shows the
# allele-balanced state sits at log2 +0.2..+0.3, the "losses" are too shallow for
# 90% purity, and the TP53/RB1 "homozygous" deletions are only log2 -1.6. FACETS
# fits purity + ploidy jointly from depth ratio AND germline-het allele fractions,
# giving integer total/minor copy number per segment, which everything downstream
# (WGD, LOH fraction, HRD score, CCF of somatic calls, cohort comparison) needs.
#
# Inputs : *.dedup.recal.bam (tumor ISM556046-2, normal ISM556054-2), hg19 with 'chr'
# Refs   : dbSNP common (GRCh37, downloaded), refGene (exon targets), VEP FASTA
# Output : results/04_copynumber/
#   facets_pileup.csv.gz         snp-pileup counts (tumor+normal) at common SNPs in exons
#   facets_cval{150,300}_*       fit per cval: segments.tsv, fit.tsv (purity/ploidy/dipLogR), plots
#   exon_depth_TP53_RB1.tsv      per-exon tumor/normal depth ratio -> which exons are deleted
#   normal_exon_depth_panel.tsv  normal-only per-exon depth on the predisposition panel (germline CNV screen)
# =============================================================================
set -euo pipefail
PROJ="${PROJ:-$PWD}"
DATA="$PROJ/data/sema4"; REF="$PROJ/refs"; OUT="$PROJ/results/04_copynumber"
mkdir -p "$REF" "$OUT"
FACETS_ENV="${FACETS_ENV:-/home/ubuntu/miniforge3/envs/facets}"
FASTA="${FASTA:-$(ls "$REF"/vep/homo_sapiens/*/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz | head -1)}"
DBSNP_URL="${DBSNP_URL:-https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz}"
REFGENE_URL="${REFGENE_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz}"
NT=$(nproc)
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

TBAM=$(find -L "$DATA" -name 'ISM556046-2.dedup.recal.bam' | head -1)
NBAM=$(find -L "$DATA" -name 'ISM556054-2.dedup.recal.bam' | head -1)
[[ -f "$TBAM" && -f "$NBAM" ]] || { log "ERROR: tumor/normal BAM not found under $DATA"; exit 1; }
log "tumor : $TBAM"; log "normal: $NBAM"
CHR_PREFIX=$(samtools view -H "$TBAM" | awk '$1=="@SQ"{sub("SN:","",$2); print ($2 ~ /^chr/) ? "chr" : ""; exit}')
log "BAM chromosome prefix: '${CHR_PREFIX}'"

# ---------- exon targets (all refGene coding exons +/-100 bp) : restricts the pileup to captured territory ----------
[[ -s "$REF/refGene_hg19.txt.gz" ]] || curl -sSL -o "$REF/refGene_hg19.txt.gz" "$REFGENE_URL"
if [[ ! -s "$REF/exons_hg19.bed" ]]; then
  zcat < "$REF/refGene_hg19.txt.gz" | awk -v OFS='\t' '$3 ~ /^chr([0-9]+|X|Y)$/ {
      n=split($10,s,","); split($11,e,",");
      for(i=1;i<n;i++){ st=s[i]-100; if(st<0)st=0; print $3, st, e[i]+100, $13 } }' \
    | sort -k1,1V -k2,2n | awk -v OFS='\t' '   # merge overlapping intervals (no bedtools on the box)
        NR==1{c=$1;s=$2;e=$3;g=$4;next}
        $1==c && $2<=e {if($3>e)e=$3; if(index(","g",", ","$4",")==0)g=g","$4; next}
        {print c,s,e,g; c=$1;s=$2;e=$3;g=$4} END{print c,s,e,g}' > "$REF/exons_hg19.bed"
fi
log "exon targets: $(wc -l < "$REF/exons_hg19.bed") intervals"

# ---------- dbSNP common SNPs, restricted to exon targets, chrom names matched to the BAM ----------
if [[ ! -s "$REF/dbsnp_common_exons.vcf.gz" ]]; then
  [[ -s "$REF/dbsnp_common_GRCh37.vcf.gz" ]] || { log "downloading dbSNP common (GRCh37, ~1.2 GB)"; curl -sSL -o "$REF/dbsnp_common_GRCh37.vcf.gz" "$DBSNP_URL"; }
  [[ -s "$REF/dbsnp_common_GRCh37.vcf.gz.tbi" ]] || tabix -f "$REF/dbsnp_common_GRCh37.vcf.gz"
  sed 's/^chr//' "$REF/exons_hg19.bed" | cut -f1-3 > "$OUT/exons_nochr.bed"
  if [[ -n "$CHR_PREFIX" ]]; then
    for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > "$OUT/chr_add.txt"; echo "MT chrM" >> "$OUT/chr_add.txt"
    bcftools view -v snps -R "$OUT/exons_nochr.bed" "$REF/dbsnp_common_GRCh37.vcf.gz" -Ou \
      | bcftools annotate --rename-chrs "$OUT/chr_add.txt" -Oz -o "$REF/dbsnp_common_exons.vcf.gz"
  else
    bcftools view -v snps -R "$OUT/exons_nochr.bed" "$REF/dbsnp_common_GRCh37.vcf.gz" -Oz -o "$REF/dbsnp_common_exons.vcf.gz"
  fi
  tabix -f "$REF/dbsnp_common_exons.vcf.gz"
fi
log "common SNPs in exon targets: $(bcftools index -n "$REF/dbsnp_common_exons.vcf.gz")"

# ---------- snp-pileup: normal first, then tumor (FACETS convention) ----------
if [[ ! -s "$OUT/facets_pileup.csv.gz" ]]; then
  log "snp-pileup (min MAPQ 15, min BAQ 20, pseudo-SNP every 100 bp) — this is the slow part"
  "$FACETS_ENV/bin/snp-pileup" -g -q15 -Q20 -P100 -r25,0 "$REF/dbsnp_common_exons.vcf.gz" \
      "$OUT/facets_pileup.csv.gz" "$NBAM" "$TBAM"
fi
log "pileup rows: $(zcat < "$OUT/facets_pileup.csv.gz" | wc -l)"

# ---------- FACETS fit at two cvals (150 = default for exomes; 300 = coarser, fewer spurious segments) ----------
log "FACETS fit"
"$FACETS_ENV/bin/Rscript" "$PROJ/scripts/lib/facets_fit.R" "$OUT/facets_pileup.csv.gz" "$OUT" 150 300

# ---------- exon-level depth at TP53 / RB1 : which exons are actually deleted? ----------
zcat < "$REF/refGene_hg19.txt.gz" | awk -v OFS='\t' -v p="$CHR_PREFIX" '
  ($13=="TP53" && $2=="NM_000546") || ($13=="RB1" && $2=="NM_000321") {
    n=split($10,s,","); split($11,e,","); strand=$4;
    for(i=1;i<n;i++){ ex = (strand=="+") ? i : n-i; c=$3; if(p=="") sub(/^chr/,"",c); print c, s[i], e[i], $13"_exon"ex } }' \
  | sort -k1,1V -k2,2n > "$OUT/tp53_rb1_exons.bed"
# exon BED with chromosome names matching the BAM
if [[ -n "$CHR_PREFIX" ]]; then cut -f1-3 "$REF/exons_hg19.bed" > "$OUT/allexons_bam.bed"; else sed 's/^chr//' "$REF/exons_hg19.bed" | cut -f1-3 > "$OUT/allexons_bam.bed"; fi
for s in tumor normal; do
  b=$([[ $s == tumor ]] && echo "$TBAM" || echo "$NBAM")
  mosdepth -t 4 -n -Q 20 --by "$OUT/tp53_rb1_exons.bed" "$OUT/exon_$s" "$b"
  # genome-wide normaliser: mean depth over ALL exon targets (summary 'total_region' only covers the --by BED)
  [[ -s "$OUT/allexons_$s.mosdepth.summary.txt" ]] && awk '$1=="total_region" && $4>0{f=1} END{exit !f}' "$OUT/allexons_$s.mosdepth.summary.txt" || mosdepth -t 8 -n -Q 20 --by "$OUT/allexons_bam.bed" "$OUT/allexons_$s" "$b"
done
TN=$(awk '$1=="total_region"{print $4}' "$OUT/allexons_tumor.mosdepth.summary.txt"); NN=$(awk '$1=="total_region"{print $4}' "$OUT/allexons_normal.mosdepth.summary.txt")
log "mean exon-target depth: tumor $TN  normal $NN"
paste <(zcat < "$OUT/exon_tumor.regions.bed.gz") <(zcat < "$OUT/exon_normal.regions.bed.gz" | cut -f5) \
  | awk -v OFS='\t' -v tn="$TN" -v nn="$NN" 'BEGIN{print "chrom","start","end","exon","tumor_depth","normal_depth","ratio_T_over_N","log2_ratio"}
      {r=($6>0)?($5/tn)/($6/nn):"NA"; print $1,$2,$3,$4,$5,$6,r,(r=="NA"||r<=0)?"NA":log(r)/log(2)}' > "$OUT/exon_depth_TP53_RB1.tsv"
log "TP53/RB1 exon depth -> $OUT/exon_depth_TP53_RB1.tsv"

# ---------- germline CNV screen: normal-only per-exon depth across the step-01 predisposition panel ----------
# No panel-of-normals, so this is coarse: a heterozygous germline exon deletion shows ~0.5 of the sample's
# own median exon coverage; a whole-gene one shows it across every exon. Anything <0.6 is flagged for a look.
PANEL="$PROJ/results/01_germline/panel_genes.txt"
if [[ -s "$PANEL" ]]; then
  zcat < "$REF/refGene_hg19.txt.gz" | awk -v OFS='\t' -v p="$CHR_PREFIX" 'NR==FNR{g[$1]=1; next}
    ($13 in g) && $3 ~ /^chr([0-9]+|X|Y)$/ { n=split($10,s,","); split($11,e,","); c=$3; if(p=="") sub(/^chr/,"",c);
      for(i=1;i<n;i++) print c, s[i], e[i], $13 }' "$PANEL" - | sort -k1,1V -k2,2n -u > "$OUT/panel_exons.bed"
  mosdepth -t 4 -n -Q 20 --by "$OUT/panel_exons.bed" "$OUT/panel_normal" "$NBAM"
  MED=$(zcat < "$OUT/panel_normal.regions.bed.gz" | cut -f5 | sort -g | awk '{d[NR]=$1} END{print (NR%2)?d[(NR+1)/2]:(d[NR/2]+d[NR/2+1])/2}')
  zcat < "$OUT/panel_normal.regions.bed.gz" | awk -v OFS='\t' -v med="$MED" 'BEGIN{print "chrom","start","end","gene","normal_depth","rel_to_median","flag"}
      {r=$5/med; print $1,$2,$3,$4,$5,sprintf("%.2f",r),(r<0.6?"LOW":(r>1.5?"HIGH":""))}' > "$OUT/normal_exon_depth_panel.tsv"
  log "panel exons flagged LOW in normal: $(grep -c 'LOW$' "$OUT/normal_exon_depth_panel.tsv" || true) of $(( $(wc -l < "$OUT/normal_exon_depth_panel.tsv") - 1 ))"
fi

rm -f "$OUT"/*.per-base.bed.gz* "$OUT"/*.mosdepth.global.dist.txt "$OUT"/allexons_*.regions.bed.gz*
log "DONE"
