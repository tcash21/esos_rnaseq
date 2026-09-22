#!/usr/bin/env bash
# chimera_dna_check.sh — do the RNA chimeras USP39–CTNNA2 and PGAP1–DNAH7 (step 03) have matching
# DNA rearrangement reads in the tumor BAM? Discordant pairs / split reads between the two gene loci.
# Run on EC2:  bash scripts/run_ec2.sh 03b   (symlinked as scripts/03b_chimera_dna_check.sh)
set -euo pipefail
PROJ="${PROJ:-$PWD}"; OUT="$PROJ/results/03_fusions"; mkdir -p "$OUT"
T=$(find -L "$PROJ/data/sema4" -name 'ISM556046-2.dedup.recal.bam' | head -1)
N=$(find -L "$PROJ/data/sema4" -name 'ISM556054-2.dedup.recal.bam' | head -1)
# hg19 gene bodies (+/- 50 kb): USP39 2:85,833,000-85,870,000; CTNNA2 2:79,593,000-80,875,000; PGAP1 2:197,697,000-197,790,000; DNAH7 2:196,602,000-196,776,000
pair(){ # <name> <regionA> <chromB> <startB> <endB>
  local name=$1 ra=$2 cb=$3 sb=$4 eb=$5
  echo "== $name : reads in $ra whose mate maps to $cb:$sb-$eb  (tumor, then normal)"
  for b in "$T" "$N"; do
    samtools view -q 20 -F 0x40E "$b" "$ra" | awk -v cb="$cb" -v sb="$sb" -v eb="$eb" -F'\t' \
      '(($7=="=" && $3==cb) || $7==cb) && $8>=sb && $8<=eb {print int($4/1000)*1000, "->", int($8/1000)*1000}' | sort | uniq -c | sort -rn | head -6
    echo "   -- total: $(samtools view -q 20 -F 0x40E "$b" "$ra" | awk -v cb="$cb" -v sb="$sb" -v eb="$eb" -F'\t' '(($7=="=" && $3==cb) || $7==cb) && $8>=sb && $8<=eb' | wc -l)"
  done
  echo "== $name : split reads (SA tag) in $ra pointing into $cb:$sb-$eb (tumor)"
  samtools view -q 20 "$T" "$ra" | awk -v cb="$cb" -v sb="$sb" -v eb="$eb" -F'\t' '/SA:Z:/{for(i=12;i<=NF;i++) if($i~/^SA:Z:/){split($i,a,":"); split(a[3],b,","); if(b[1]==cb && b[2]>=sb && b[2]<=eb) print $4, $6, "->", b[2], b[3]}}' | sort | uniq -c | sort -rn | head -8
}
{
pair "USP39-CTNNA2" 2:85783000-85920000 2 79543000 80925000
pair "CTNNA2-USP39" 2:79543000-80925000 2 85783000 85920000
pair "PGAP1-DNAH7"  2:197647000-197840000 2 196552000 196826000
pair "DNAH7-PGAP1"  2:196552000-196826000 2 197647000 197840000
# MAP4K4 3' region (kataegis cluster at hg19 2:102,476,129-102,477,440) x the RNA junction acceptor seen by STAR at
# GRCh38 2:102,453,450 (~hg19 2:103.07 Mb). 173 RNA split reads shared one acceptor -> expect a DNA breakpoint here.
pair "MAP4K4-3prime x 103.07Mb" 2:102420000-102530000 2 102950000 103200000
pair "103.07Mb x MAP4K4-3prime" 2:102950000-103200000 2 102420000 102530000
} | tee "$OUT/chimera_dna_check.txt"
