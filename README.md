# mom_esos — extraskeletal osteosarcoma (ESOS) case genomics

Tumor–normal whole-exome + tumor whole-transcriptome (Sema4 Signal WES/WTS, 2020, hg19).
Tumor purity 90%, tumor 276x, normal 122x, RNA 327M reads. Diagnosis: high-grade ESOS, abdominal mass.
Reported drivers: homozygous deletions of TP53 (17p13.1) and RB1 (13q14.2); TMB 3.77; MSS; no fusions reported.

## Layout
```
scripts/NN_*.sh, NN_*.R   analysis steps, numbered; each writes results/NN_*/
scripts/run_ec2.sh        runs a step on the EC2 instance and pulls results back
infra/                    AWS setup, bootstrap and Dropbox->EC2 transfer scripts
data/sema4/               Sema4 delivery (NOT in git)
refs/                     downloaded references (NOT in git)
results/                  step outputs + logs (NOT in git)
```

## Running a step
Compute happens on the EC2 instance (`i-09847e2ec5b85a3bc`, bucket `esos-sema4-576d839b`):
it already holds the data, has bcftools, and can reach ClinVar/gnomAD/UCSC.
```bash
bash scripts/run_ec2.sh 01        # ~10-20 min; results land in results/01_germline/
cat results/01_germline/germline_report.md
```
`esos-stop` (zsh helper) when done — the disk is billed while the instance runs.

## Roadmap
| Step | Question | Inputs (already in delivery) | Status |
|---|---|---|---|
| 01 germline re-check | Any inherited risk missed in 2020? | `*.germline.vcf` (normal) | scripts ready |
| 02 expression | Does the tumor look like bone OS or soft-tissue sarcoma? What is over/under-expressed (KDM3A, CDK/cyclins, immune)? | `*.genes.results` (RSEM) vs TARGET-OS + TCGA-SARC (UCSC Xena) | next |
| 03 fusions | Anything in the full candidate lists the report dropped? | STAR-Fusion / FusionCatcher / FusionInspector TSVs | |
| 04 copy number | Ploidy, whole-genome doubling, LOH, HRD scar score; confirm TP53/RB1 homozygous loss | GATK `.seg` now; FACETS on BAMs later (EC2) | |
| 05 somatic | Re-annotate Mutect2 calls vs current databases; TP53/RB1 vs OS cohorts (GENIE) | `*.mutect2.twicefiltered.vcf.gz` | |
| 06 immune / HLA | Immune infiltrate (quanTIseq/MCP-counter), HLA type, HLA LOH | RNA `genes.results`; BAMs (EC2) | |
| 07 deposit + write-up | Zenodo (processed, de-identified) + Molecular Case Studies submission | everything above | |

## Data sources (all free unless noted)
- ClinVar, gnomAD v2.1.1 (GRCh37) — variant interpretation
- UCSC Xena: TARGET-OS and TCGA-SARC RSEM expression — comparison cohorts
- cBioPortal / AACR GENIE — osteosarcoma mutation & CNA frequencies
- CIViC (open), OncoKB (academic license free; commercial use requires license)
- COSMIC signatures via SigProfiler (free for non-commercial)
- Paid, optional, high value: methylation array (EPIC) on the archived tumor block → DKFZ sarcoma classifier (~$300–600)

## Rules
- Never commit `data/`, `refs/`, `results/`, or anything derived from the germline VCF.
- Everything is research-grade. Germline findings go to a genetic counselor before anyone acts on them.
