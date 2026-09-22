#!/usr/bin/env bash
# make_review_bundle.sh — assemble a de-identified bundle for external reviewers.
#   bash scripts/make_review_bundle.sh            -> review_bundle_<date>.tar.gz in the project root
# Includes: scripts/, README, NOTEBOOK, REVIEW, manuscript/, results/02-06 reports/tables/plots, logs.
# Excludes: data/ (delivery incl. the Sema4 PDF with identifiers), refs/, results/01_germline/ (germline),
#           results/06_hla/hla_* (HLA alleles), BAM/VCF/pileup/FASTQ-derived bulk files, the 890 MB chimeric junction file.
set -euo pipefail
cd "$(dirname "$0")/.."
stamp=$(date +%Y%m%d); out="review_bundle_${stamp}.tar.gz"
tar -czf "$out" \
  --exclude='results/01_germline' --exclude='results/06_hla/hla_*' --exclude='results/06_hla/opti_*' \
  --exclude='*.vcf' --exclude='*.vcf.gz' --exclude='*.vcf.gz.tbi' --exclude='*.bam' --exclude='*.bai' \
  --exclude='*pileup*' --exclude='*.Chimeric.out.junction' --exclude='*.regions.bed.gz*' --exclude='*.per-base*' \
  --exclude='*.ReadsPerGene.out.tab' --exclude='*.SJ.out.tab' --exclude='*_STARgenome' --exclude='*_STARpass1' \
  --exclude='.git' --exclude='.DS_Store' \
  README.md NOTEBOOK.md REVIEW.md manuscript scripts results
echo "wrote $out ($(du -h "$out" | cut -f1))"
echo "contents check — these must be ABSENT:"; tar -tzf "$out" | grep -E '01_germline|hla_|Final Lab|\.vcf|\.bam' || echo "  none found (good)"
