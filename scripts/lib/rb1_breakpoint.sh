#!/usr/bin/env bash
# rb1_breakpoint.sh — look for split (SA-tag) and discordant reads at the RB1 exon 18-24 deletion
# boundaries in the tumor BAM. Writes results/04_copynumber/rb1_breakpoint.txt. Run on EC2:
#   bash scripts/run_ec2.sh 04b   (this file is symlinked as scripts/04b_rb1_breakpoint.sh)
set -euo pipefail
PROJ="${PROJ:-$PWD}"; OUT="$PROJ/results/04_copynumber"; mkdir -p "$OUT"
T=$(find -L "$PROJ/data/sema4" -name 'ISM556046-2.dedup.recal.bam' | head -1)
{
echo "== split reads (SA tag), MAPQ>=20, in 13:48950000-49000000 (exon 17 -> exon 18 region); primary pos, CIGAR -> supplementary pos"
samtools view -q 20 "$T" 13:48950000-49000000 \
  | awk -F'\t' '/SA:Z:/{for(i=12;i<=NF;i++) if($i~/^SA:Z:/){split($i,a,":"); split(a[3],b,","); print $3":"$4, $6, "->", b[1]":"b[2], b[3]}}' \
  | sort | uniq -c | sort -rn | head -15
echo
echo "== discordant same-chrom pairs (mate >20 kb away), reads starting in 13:48980000-48992000, binned to 1 kb"
samtools view -q 20 -F 0x40E "$T" 13:48980000-48992000 \
  | awk -F'\t' '$7=="=" && ($8-$4>20000 || $4-$8>20000){print int($4/1000)*1000, "->", int($8/1000)*1000}' | sort | uniq -c | sort -rn | head -10
echo
echo "== same, reads starting in 13:49050000-49062000 (exon 24/25 side)"
samtools view -q 20 -F 0x40E "$T" 13:49050000-49062000 \
  | awk -F'\t' '$7=="=" && ($8-$4>20000 || $4-$8>20000){print int($4/1000)*1000, "->", int($8/1000)*1000}' | sort | uniq -c | sort -rn | head -10
echo
echo "== soft-clipped reads (>=20 bp clip), MAPQ>=20, 13:48985000-48988000 by position"
samtools view -q 20 "$T" 13:48985000-48988000 | awk -F'\t' '$6 ~ /[0-9]{2,}S/ {print $4, $6}' | sort | uniq -c | sort -k2,2n | awk '$1>=2' | head -20
} | tee "$OUT/rb1_breakpoint.txt"
