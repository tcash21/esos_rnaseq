# Review request — ESOS case genomics (mom_esos)

Prepared 2026-09-22 for external review by a bioinformatics team. Everything here is research-grade.
Please read this file, then `NOTEBOOK.md` (findings with confidence tags, corrections, run log), then the
per-step reports listed in §3. Scripts are in `scripts/`; each step's header explains what it does and why.

## 1. The case and the question

Tumor–normal whole-exome (tumor 276×, normal 122×) plus tumor total-RNA sequencing (327 M reads) of a
high-grade **extraskeletal osteosarcoma** (ESOS) of the **retroperitoneum** in a woman in her mid-60s,
sequenced clinically in 2020 (Sema4 Signal WES/WTS, hg19). The clinical report called homozygous TP53 and
RB1 deletions, TMB 3.77, MSS, no fusions, no hereditary findings.

The research question: **is this ESOS biologically a bone osteosarcoma arising in soft tissue, or a
soft-tissue sarcoma with osteoid?** That distinction decides which chemotherapy family is used and has
never been settled by a trial (ESOS is ~1% of soft-tissue sarcomas). Secondary aims: what did the
clinical report miss or not attempt, and which findings are worth depositing and publishing.

## 2. Headline claims we are asking you to test

| # | Claim | Our confidence | Where |
|---|---|---|---|
| A | Purity 0.57, ploidy 3.4, whole-genome doubled (FACETS; Sema4 assumed 90% purity, diploid) | solid (one caller, two cvals agree) — **wants a 2nd caller** | `results/04_copynumber/` |
| B | TP53: whole-gene homozygous deletion. RB1: somatic p.Thr502Ile on both retained copies + copy-neutral LOH of exons 1–17 + homozygous deletion from intron 17 (junction chr13:48,986,411 ↔ 49,237,832, inverted, 36 split / 94 discordant reads) | solid | `04`, `05`, `rb1_breakpoint.txt` |
| C | Transcriptome clusters with bone OS (TARGET-OS, n=88), not with any TCGA-SARC subtype (n=259), on a like-for-like STAR/GENCODE v36 quantification | solid | `results/02_expression/gdc_requant/` |
| D | 18-copy 2p11.2 amplicon (~5 Mb, hg19 chr2:83.08–88.13 Mb) containing KDM3A, transcribed across its span; two RNA chimeras coincide with its DNA boundaries | likely | `04`, `02 gdc`, `03` |
| E | MAP4K4: four clonal strand-coordinated C>T/C>G SNVs in 1.3 kb (3/4 TpC), on ~3 of 5 copies (pre-WGD), adjacent (~10–15 kb) to an RNA-evidenced breakpoint (MAP4K4–SLC9A2, 173 STAR split reads at one acceptor). Interpreted as APOBEC kataegis at a rearrangement, not four drivers. Sema4 reported one as a VUS | likely | `05`, `03` |
| F | HRD scar score 75–81 (Telli definitions on exome segments) with clean HR genes germline and somatic | weak — **method needs your check** | `04` |
| G | Microenvironment: myeloid/macrophage-dominated (MCP-counter monocytic score above all 88 TARGET-OS), CD8 at OS median, PD-L1 low; HLA class I heterozygous, all alleles expressed, no HLA LOH | likely | `06` |
| H | No limb HOX code (HOXA9/10/11, HOXC10, PITX1 near zero vs TARGET-OS). Retroperitoneal site is lumbar-level where HOX10 paralogs are normally expressed, so this may reflect cell of origin | weak — needs the right control | `02 gdc` |
| I | Germline: no P/LP variant in 51 predisposition genes (ACMG + sarcoma genes), no germline exon/gene deletion in TP53/RB1 (coarse screen); one incidental carrier-type variant outside the panel | solid / coarse | local only, see §5 |
| J | Sema4's somatic VCF `PASS` set (1,177) is ~90% Mutect2-rescued noise (<5% AF read clusters); 87 calls survive an explicit policy; TMB 1.8/Mb by that policy | solid | `05` |

## 3. Specific items to check (numbered as in NOTEBOOK §4)

1. **FACETS solution** — purity 0.57 / ploidy 3.42 / dipLogR −0.49, cval 150 and 300. FACETS flagged and rejected a near-diploid alternative ("mafR larger than expected if −0.34 is diploid"). Please eyeball `facets_cval150.png` and, if you can, run Sequenza or ASCAT on the same pileup/BAMs (`04_copynumber.sh` documents the snp-pileup call). Everything downstream (CCF, multiplicity, WGD, HRD, timing) inherits this.
2. **HRD score implementation** — `scripts/04_copynumber.R` re-implements LOH/TAI/LST (Telli 2016) on exome segments; not scarHRD. Are the TAI and LST definitions right? Is 75–81 credible for an exome after WGD?
3. **RB1 junction** — `results/04_copynumber/rb1_breakpoint.txt`. Inverted intrachromosomal junction: fold-back/BFB vs. a two-breakpoint complex event? Exons 25–27 read ~1 copy on per-exon depth but FACETS says 0 through 49.28 Mb.
4. **Germline artifact triage** — four ClinVar "pathogenic" hits dismissed on alt-fraction / mapping-quality grounds (PRSS1 ×2 at AF 0.20, SLC9B1 at 0.17, GSTT2 at 5 reads). Agree? (Local report, §5.)
5. **Germline CNV screen** — normal-BAM exon depth vs sample median, no panel of normals. POT1 came out at 0.64× with patchy exons; we read it as capture variability. Sensitivity for single-exon events is poor — is there a better approach without a PoN?
6. **Somatic filter policy** — `scripts/05_somatic.R` header. PASS, ≥8 alt reads, AF ≥0.05 (SNV) / ≥0.10 (indel), normal AF ≤0.02; Sema4's `lowAfT`/`readsInN` flags replaced by these thresholds; `mutectFiltOverride` deliberately *not* excluded (Mutect2's clustered-events filter removes kataegis). Reasonable?
7. **RB1 p.Thr502Ile** (c.1505C>T; COSMIC/HGMD ids) — unreported by Sema4 at 50% AF, 117×. Classification, and whether it belongs on a corrected clinical record.
8. **MAP4K4 cluster = kataegis?** Four events is small for a formal call. Any MAP4K4 sarcoma literature?
9. **TMB** — 1.8/Mb (our policy) vs 3.77 (Sema4, definition unknown). Only wording matters.
10. **Expression pipeline** — tumor re-quantified with the GDC recipe (STAR 2-pass, GRCh38, GENCODE v36, unstranded counts → TPM; `02b_requant_rna.sh`) to match Xena's GDC-hub cohorts. Remaining confounder: tumor is stranded **total RNA** (51% intronic bases), cohorts are poly-A. We excluded histone/sn/scaRNA/7SL/MT and PAR1 genes from genome-wide lists. Is the rank-based comparison (quantile mapping + Spearman + centroid + PCA projection) adequate, and are single-gene z-scores defensible with a library-prep mismatch?
11. **Breakpoint chimeras** — USP39–CTNNA2, PGAP1–DNAH7, KDM3A–CYTOR, MAP4K4–SLC9A2 as transcribed DNA rearrangements. Exome DNA could not confirm (partner sites outside capture). Worth WGS or junction PCR?
11b. **HOX / positional identity** — see claim H. The right control is HOX10 expression in TCGA-SARC retroperitoneal tumors (DDLPS/LMS) vs limb tumors; not yet done.
12. **Sema4 CN files** — the delivery has `segData.seg`, `.female.seg`, `.male.seg`; we used the unsuffixed one for comparison. Confirm it is the reported one.
13. **PPP1CB–ALK** — 43 STAR chimeric reads; dismissed (ALK 0.4 TPM, six scattered breakpoints, non-canonical). Please confirm, since ALK is actionable.

## 4. Known errors, fixed, and limitations (so you don't re-find them)

See NOTEBOOK §3 (corrections table) and §5. In short: runner/SSM quirks; an idle-stop alarm that killed a download; Sema4 VCF header lines bcftools rejects; GT written as bare `1` with zygosity in `INFO/CALL_TYPE`; `chr` in VCFs but not BAMs; a mosdepth normaliser bug in the first exon-depth run; HMGA2/VEGFA/BGLAP/CDK4 z-scores withdrawn after re-quantification; "T-cell-poor" corrected to "myeloid-dominated, CD8 at median" after deconvolution.

## 5. What you will and will not receive (PHI / consent)

The repository (scripts, NOTEBOOK, REVIEW, manuscript draft) contains no patient identifiers and no germline
variant detail. The reviewer bundle (`scripts/make_review_bundle.sh`) adds `results/02–06` tables, plots and
reports and the somatic VCF-derived tables. It **excludes**: the Sema4 PDF (name, DOB, MRN), everything in
`results/01_germline/` (germline variants), the HLA allele table, BAMs/FASTQs, and the raw delivery. If the
review needs germline-level detail (items 4, 5), that will be shared separately under the family's consent.

## 6. How to reproduce

- EC2 steps (need BAMs): `bash scripts/run_ec2.sh <00|01|02b|03b|04|04b|05|06>`; instance/bucket in `scripts/run_ec2.sh`.
- Mac steps (need only tables + Xena downloads, which the scripts fetch into `refs/xena/`):
  `PROJ=$PWD QUANT=gdc Rscript scripts/02_expression.R`, then `03_fusions.R`, `04_copynumber.R`, `05_somatic.R`, `06_immune.R`.
- Software: bcftools/htslib 1.21–1.24, VEP 116 (GRCh37 cache), FACETS 0.6.2 + snp-pileup, mosdepth, STAR 2.7.10b, OptiType 1.5, R 4.4 with data.table, ggplot2, MCPcounter 1.2, quantiseqr 1.14. dbSNP common b151, ClinVar 2026-09-13, gnomAD v2.1.1, GENCODE v36, Xena GDC hub (TARGET-OS, TCGA-SARC star_tpm).

## 7. What would most change the write-up

1. A second purity/ploidy caller disagreeing with FACETS.
2. The HRD score collapsing under a proper implementation.
3. A reason the OS-vs-STS transcriptome call could be a library-prep artifact (we don't see one: an artifact should push the sample away from both cohorts, not toward one).
4. Any real fusion we scored as an artifact.
