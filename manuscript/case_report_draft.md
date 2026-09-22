# Integrated genomic and transcriptomic profiling of a retroperitoneal extraskeletal osteosarcoma reveals a bone-osteosarcoma genome and transcriptome in a soft-tissue site

*Draft v0.1 — 2026-09-22 — for internal and expert review. Research-grade; not for clinical use.*
*Target: Cold Spring Harbor Molecular Case Studies (Research Report). Placeholders in [brackets] need clinical records or consent. References marked [verify] were cited from memory and must be checked.*

Authors: [Tanya Cashorali]¹, [reviewing bioinformatics team], [treating clinician(s), if participating], …
¹ TCB Analytics, [city]

Corresponding author: [ ]

---

## Abstract

Extraskeletal osteosarcoma (ESOS) is a rare, high-grade sarcoma that produces osteoid in soft tissue without skeletal attachment. It is classified among soft-tissue sarcomas yet is histologically indistinguishable from conventional osteosarcoma of bone, and whether it should be treated with bone-osteosarcoma or soft-tissue-sarcoma chemotherapy regimens remains unresolved; no randomized trial exists and molecular characterization at case level is sparse. We report integrated tumor–normal whole-exome and tumor whole-transcriptome analysis of a retroperitoneal ESOS in a woman in her mid-60s, re-analyzing a 2020 clinical dataset with current tools and cohort comparisons. Allele-specific copy-number analysis showed a whole-genome-doubled tumor (ploidy 3.4, purity 0.57; the clinical estimate was 90%) with biallelic inactivation of TP53 (whole-gene homozygous deletion) and RB1 by three events (a clonal p.Thr502Ile missense on both retained copies, copy-neutral loss of heterozygosity, and a homozygous intragenic deletion from intron 17 with an inverted junction resolved at base-pair resolution), together with 17p11.2 and 8q gain and a novel 18-copy 2p11.2 amplicon encompassing KDM3A that was transcriptionally active across its span. The mutational landscape was quiet (1.8 non-synonymous mutations/Mb) except for a clonal APOBEC-type kataegis cluster in MAP4K4 adjacent to an RNA-evidenced rearrangement breakpoint, which the clinical report had partly captured as a single variant of uncertain significance. No gene fusion was present. Against 88 TARGET-OS bone osteosarcomas and 259 TCGA soft-tissue sarcomas quantified with the same pipeline, the tumor's transcriptome clustered unambiguously with bone osteosarcoma (median Spearman ρ 0.64 vs ≤0.42 for any soft-tissue subtype; all 25 nearest neighbours were osteosarcomas) while lacking the limb HOX code of skeletal osteosarcoma. The microenvironment was macrophage-dominated with median CD8 T-cell content and intact HLA class I. Germline analysis found no cancer-predisposition variant. This case supports the view that ESOS is osteosarcoma by genome and transcriptome, not merely by histology, and illustrates what re-analysis of clinical sequencing data can add for rare tumors.

**Keywords:** extraskeletal osteosarcoma; osteosarcoma; TP53; RB1; whole-genome doubling; KDM3A; kataegis; RNA-seq; TARGET

---

## Introduction

Extraskeletal osteosarcoma accounts for roughly 1% of soft-tissue sarcomas and 2–4% of all osteosarcomas [verify: Longhi 2017; WHO 2020]. It arises most often in the thigh and retroperitoneum, in patients decades older than those with conventional osteosarcoma of bone, and carries a worse prognosis [verify]. Because it is defined by osteoid production in a soft-tissue location, its diagnosis is histological; but its biological identity — an osteosarcoma that happens to arise outside bone, or a soft-tissue sarcoma with osteoblastic differentiation — has direct consequences for management, since bone osteosarcoma is treated with methotrexate/doxorubicin/cisplatin-based regimens and soft-tissue sarcomas with doxorubicin/ifosfamide. Retrospective series disagree on whether osteosarcoma-type chemotherapy benefits ESOS [verify: Longhi 2017; Paludo 2018], and the disease is too rare for a randomized comparison.

Conventional osteosarcoma has a characteristic genome: near-universal inactivation of TP53 and frequent loss of RB1, extensive copy-number instability with chromothripsis and whole-genome doubling, recurrent gains at 6p12–21, 8q24, 17p11.2 and 19q12, kataegis, and a striking absence of recurrent driver point mutations or fusions [Chen 2014; Perry 2014; Behjati 2017; Kovac 2015]. Whether ESOS shares this landscape is known only from small series [verify], and to our knowledge no ESOS case has been reported with paired exome, deep transcriptome and cohort-level transcriptomic comparison.

Here we re-analyzed a clinical tumor–normal exome and transcriptome dataset generated in 2020 for a patient with retroperitoneal ESOS, using current annotation resources, allele-specific copy-number modelling, and like-for-like comparison with public osteosarcoma and soft-tissue-sarcoma cohorts. We asked what the clinical report had established, what it had not attempted, and whether the tumor's molecular identity is that of bone osteosarcoma.

---

## Results

### Clinical presentation

A woman in her mid-60s presented with [symptoms; date] and was found to have a [size] retroperitoneal mass involving [organs/vessels]. Biopsy showed a high-grade sarcoma with malignant osteoid production and no attachment to bone on imaging, diagnosed as high-grade extraskeletal osteosarcoma [pathology details; IHC if available]. [Staging: metastatic disease at diagnosis yes/no.] The tumor was judged unresectable and she was treated with radiation alone [dose/fractions/dates]; chemotherapy was [not recommended / declined] [reason, if documented]. [Follow-up and outcome.] Clinical tumor–normal whole-exome and tumor whole-transcriptome sequencing (Sema4 Signal WES/WTS) was performed on the diagnostic tissue block and peripheral blood in mid-2020; the report identified homozygous deletions of TP53 and RB1, a tumor mutational burden of 3.77 mutations/Mb, microsatellite stability, one variant of uncertain significance (MAP4K4 p.E589Q), no fusions, and no hereditary findings. The family later obtained the underlying data for research re-analysis [consent statement].

### Sequencing data

Tumor and normal exomes reached mean target coverages of 276× and 122×; tumor RNA (stranded total-RNA library) yielded 327 million reads. Pathologist-estimated tumor cellularity was 90%.

### Germline

Re-annotation of the normal-sample calls against ClinVar (2026-09-13), gnomAD v2.1.1 and Ensembl VEP 116 found no pathogenic or likely pathogenic variant in any of 51 cancer-predisposition genes (ACMG secondary-findings cancer genes plus RECQL4, BLM, WRN, CDKN2A, POT1 and other sarcoma-associated genes), no rare protein-altering variant of uncertain significance in these genes, and no evidence of a heterozygous germline deletion of TP53 or RB1 on exon-level read depth. Four ClinVar-pathogenic calls outside the panel were sequencing artifacts (allele fraction 0.17–0.20 at high depth in paralog-rich loci, or ≤5 reads). [One incidental heterozygous carrier-type variant outside the predisposition panel was identified and communicated for genetic counselling; details are withheld here pending consent.] The tumor's TP53 and RB1 losses are therefore somatic.

### Purity, ploidy and allele-specific copy number

Allele-specific copy-number analysis with FACETS on 1.75 million common-SNP positions gave a tumor purity of 0.57 and ploidy of 3.4, stable across segmentation stringencies (Figure 1A; Supplementary Table S1). More than half of the autosomal genome carried a major allele at ≥2 copies, meeting the criterion for whole-genome doubling [Bielski 2018], and 32% of the genome showed loss of heterozygosity. The clinical pipeline had assumed a diploid baseline and 90% purity; its reported copy numbers (e.g., "×1.5 losses" and "×2.5–3 gains") correspond, after re-fitting, to two-copy and ploidy-level states respectively. Focal amplification of 2p11.2 (hg19 chr2:83.08–88.13 Mb) reached 18 copies and contained KDM3A (Figure 1B); 17p11.2 (AURKB, COPS3) was at 7 copies with LOH, and 8q (MYC) at 6 copies. RUNX2, VEGFA, CCND3 and CCNE1 were at ploidy level (4 copies) rather than gained. Single-copy loss with LOH affected PTEN, TSC2 and CDKN1C; copy-neutral LOH affected ATRX, STK11 and CHD5. An exome-approximated homologous-recombination-deficiency scar score (LOH 19–20, TAI 31–33, LST 25–28; sum 75–81) exceeded the 42-point threshold used in breast and ovarian cancer [Telli 2016], although whole-genome doubling inflates such scores and the estimate is provisional.

### TP53 and RB1

TP53 was homozygously deleted across all 11 exons within a 46-kb focal deletion inside a 17p LOH block; the tumor:normal log-ratio (−2.1) matched the expectation for zero tumor copies at 57% purity (Figure 2A). RB1 was inactivated by three events (Figure 2B). A somatic missense variant, c.1505C>T p.Thr502Ile in the pocket A domain (COSMIC and HGMD entries [verify ids]), was present at allele fraction 0.50 (58/117 reads; 0/55 in the normal), corresponding to both retained copies; exons 1–17 were retained at two copies with loss of heterozygosity; and exons 18 onward were homozygously deleted. Split and discordant reads resolved the deletion junction to chr13:48,986,411 in intron 17, joined in inverted orientation to chr13:49,237,832, ~180 kb beyond the 3' end of RB1 (36 split reads, 94 discordant pairs), indicating a fold-back-type rearrangement rather than a simple deletion. Multiplicity analysis placed the RB1 point mutation and the TP53 and RB1 deletions before whole-genome doubling. Residual TP53 and RB1 RNA (4.4 and 21 TPM) was attributable to the 43% non-tumor cell fraction. The clinical report had correctly called both deletions but did not report the RB1 missense variant.

### Somatic point mutations and a kataegis cluster

The clinical somatic call set (1,177 PASS calls) was dominated by low-allele-fraction calls that Mutect2 had filtered and the vendor pipeline had rescued; after an explicit filter (≥8 alternate reads, allele fraction ≥0.05 for SNVs and ≥0.10 for indels, ≤0.02 in the normal), 87 somatic calls remained (82 SNVs, 5 indels), corresponding to 1.8 non-synonymous mutations/Mb. Fifty-four were clonal, and 28 had multiplicity ≥2 in gained regions, i.e. predated genome doubling. The substitution spectrum (C>T 42, C>G 14, C>A 9 of 82) was consistent with age- and APOBEC-related processes but too small for formal signature decomposition. No somatic mutation affected BRCA1, BRCA2, PALB2, RAD51C/D, ATM, ATRX or other homologous-recombination or telomere-maintenance genes.

The only clustered event was in MAP4K4: four missense SNVs within 1.3 kb (p.Glu503Lys, p.Ala507Thr, p.Glu517Gln, p.Glu620Gln; the clinically reported p.Glu589Gln corresponds to one of these under a different transcript), all G>A or G>C on the plus strand (i.e., strand-coordinated C>T/C>G), three of four in a TpC context, each at allele fraction 0.47–0.48 with no reads in the normal, and each on approximately three of the locus's five copies (Figure 3A). Three had been filtered by Mutect2 as clustered events. The independent STAR alignment of the RNA revealed 173 chimeric reads joining MAP4K4 exons to a single acceptor position ~600 kb downstream (within SLC9A2, which is not expressed), consistent with transcription across a genomic rearrangement breakpoint located within ~10–15 kb of the mutation cluster (Figure 3B). We interpret the cluster as APOBEC-type kataegis at a rearrangement junction [Nik-Zainal 2012], an established feature of osteosarcoma genomes [Behjati 2017], rather than four independent driver mutations.

### Fusions

Sema4's finding of no clinically significant fusion was upheld. Of 176 FusionCatcher candidates, none involved a sarcoma-defining partner (EWSR1, FUS, SS18, NR4A3, CIC, BCOR, NTRK1–3, HEY1–NCOA2) beyond artifact-level support; 43 were template-switching chimeras of COL1A1, the tumor's most abundant transcript; and the single FusionInspector-validated call (SEC31A–JAK2) had four junction reads, no spanning pairs and non-canonical splice sites. An independent STAR chimeric-junction analysis on 870,000 distant split reads confirmed the absence of canonical fusions and specifically excluded an ALK rearrangement (ALK 0.4 TPM; scattered non-canonical junctions). Several chimeric transcripts instead marked DNA copy-number boundaries: USP39–CTNNA2 (detected by both callers) joined the interior of the 2p11.2 amplicon to its flank, KDM3A–CYTOR was an intra-amplicon junction, and PGAP1–DNAH7 coincided with 2q segment boundaries (Figure 1B). Exome data could not confirm these junctions at the DNA level because the partner sites lie outside capture.

### Transcriptome: an osteosarcoma expression program without a limb HOX code

To compare the tumor with public cohorts on equal footing, tumor RNA was re-aligned and quantified with the GDC recipe (STAR two-pass, GRCh38, GENCODE v36; 163.7 million read pairs, 86.8% uniquely mapped) and compared with 88 TARGET-OS primary bone osteosarcomas and 259 TCGA-SARC primary soft-tissue sarcomas (leiomyosarcoma 101, dedifferentiated liposarcoma 57, undifferentiated pleomorphic sarcoma 49, myxofibrosarcoma 25, MPNST 9, synovial sarcoma 10, other 8) from the UCSC Xena GDC hub, using rank-based methods to limit residual pipeline and library-preparation effects. Over the 2,000 most variable genes, the tumor's median Spearman correlation was 0.64 with TARGET-OS versus 0.42 with the closest soft-tissue subtypes (myxofibrosarcoma, UPS) and 0.36 with leiomyosarcoma; all 25 most similar cohort samples were osteosarcomas; centroid correlation was 0.75 for osteosarcoma versus 0.52 for the next group; and in principal-component space the tumor projected inside the osteosarcoma cluster (Figure 4A). The osteoblastic program (SP7, IBSP, SPP1, ALPL, COL1A1, SATB2) was at typical osteosarcoma levels and 2–6 standard deviations above soft-tissue sarcomas; RUNX2 was at the 100th percentile of both cohorts without being amplified.

Genome-wide, the most over-expressed genes relative to osteosarcoma included six neighbours of KDM3A within the 2p11.2 amplicon (TGOLN2, ELMOD3, GGCX, MAT2A, POLR1A, PTCD3), indicating coordinate transcription across the amplified segment; KDM3A itself was 4.3 SD above the osteosarcoma mean (Figure 4B). Conversely, the posterior HOX genes that pattern the limbs and lumbosacral axis (HOXA9, HOXA10, HOXA11, HOXC10, PITX1) and MSX1 were 2.5–4.6 SD below osteosarcoma and near zero in absolute terms (Figure 4C). Smooth-muscle genes (MYOCD, ACTA2) were elevated relative to osteosarcoma but typical of soft-tissue sarcoma, consistent with retroperitoneal stroma in the non-tumor fraction. TERT was not expressed; ATRX was expressed above osteosarcoma levels while DAXX was reduced. PDGFRA and FGFR1 were at the 98th–99th percentile of osteosarcoma.

### Immune microenvironment and HLA

Deconvolution with MCP-counter and quanTIseq, run jointly on the tumor and both cohorts, showed a myeloid-dominated microenvironment: the monocytic-lineage score exceeded that of every TARGET-OS sample, neutrophil and NK-cell scores were at the 91st–94th percentiles, and quanTIseq estimated macrophages at ~17% of cells (M1 6%, M2 11%; osteosarcoma medians 1% and 7%). CD8 T-cell content was at the osteosarcoma median, CD4-type T cells were present, and CD274 (PD-L1) was low (Figure 5). HLA class I typing from normal DNA, tumor DNA and tumor RNA gave identical heterozygous genotypes at HLA-A, -B and -C, all six alleles were expressed, and the HLA locus retained both parental alleles (3 total / 1 minor copy), so antigen presentation was intact.

---

## Discussion

Every layer of this tumor's molecular profile is that of conventional osteosarcoma. The genome shows the canonical osteosarcoma configuration — biallelic TP53 and RB1 inactivation by structural events, whole-genome doubling, extensive LOH, 17p11.2 and 8q gain, kataegis at a rearrangement, and no fusion or recurrent driver point mutation [Chen 2014; Behjati 2017] — and the transcriptome clusters with bone osteosarcoma rather than with any soft-tissue sarcoma subtype when compared on a common pipeline. The dedifferentiated liposarcoma alternative, the main mimic for a retroperitoneal bone-forming sarcoma, is excluded by the absence of MDM2/CDK4 amplification and by expression. To the extent one case can speak to a classification debate, it argues that ESOS is osteosarcoma by identity and not only by histology, and it adds a molecularly characterized case to a literature in which ESOS genomics remains sparse [verify].

Three features are case-specific and, we think, of interest beyond this patient. First, an 18-copy amplicon at 2p11.2 — not among the recurrent osteosarcoma amplicons — was transcribed across its ~5-Mb span, with KDM3A, an H3K9 demethylase implicated in hypoxia response and in Ewing sarcoma biology [verify: Parrish 2015], the most compelling candidate target; chimeric transcripts at its boundaries and interior corroborate a rearranged amplicon structure. Second, the MAP4K4 cluster illustrates a practical pitfall: a clinical pipeline reported one of four adjacent clonal mutations as a variant of uncertain significance, whereas their strand coordination, TpC context, identical allele fractions and proximity to a rearrangement breakpoint identify them as a single APOBEC kataegis event with no implication of MAP4K4 as a driver. Third, the tumor lacks the limb HOX code that skeletal osteosarcomas carry as a memory of their site of origin. Because the retroperitoneum is at lumbar level, where HOX10 paralogs are normally expressed in axial mesenchyme, this absence is not explained by anatomical position alone and may reflect the tumor's cell of origin — a testable hypothesis with the retroperitoneal tumors in TCGA-SARC as a control.

The high HRD scar score without a homologous-recombination gene lesion echoes the "BRCAness" phenotype described in about a third of osteosarcomas [Kovac 2015]. We regard it as provisional: the score was derived from exome segments and whole-genome doubling inflates its components. Together with the age-associated substitution spectrum and the absence of a telomerase transcript with retained ATRX, these observations delineate the tumor's mutational processes without identifying an actionable target; the myeloid-rich, CD8-median, PD-L1-low microenvironment with intact HLA is the pattern that has not responded to checkpoint blockade in osteosarcoma trials [verify], whereas macrophage-directed and B7-H3-directed strategies would be the rational research directions.

Re-analysis also showed what a 2020 clinical report could and could not offer. The vendor called the TP53 and RB1 deletions correctly but under a diploid model with an overestimated purity, so all absolute copy numbers were wrong and copy-neutral LOH was invisible; it did not report a clonal, database-listed RB1 missense variant; its somatic VCF's PASS set was ~90% rescued noise and unusable without re-filtering; and its expression data were never analysed because the assay was not validated for expression. None of these affected the clinical decision — the tumor was unresectable and was treated with radiation — but they matter for anyone re-using clinical sequencing files for research.

Limitations are those of a single case and of exome-plus-RNA data: purity and ploidy rest on one caller pending confirmation; the RNA library was total-RNA whereas the comparison cohorts were poly-A, which we mitigated with rank-based methods and by excluding library-sensitive gene classes; rearrangement junctions inferred from RNA could not be confirmed in exome DNA; and germline copy-number screening lacked a panel of normals. [Methylation-based classification on the archived block would provide an orthogonal identity call and is planned/was not possible.]

---

## Methods

**Samples and clinical sequencing.** Tumor tissue (FFPE block, retroperitoneal mass) and peripheral blood were sequenced in 2020 by Sema4 (Signal WES/WTS, CLIA-certified): Twist exome capture, NovaSeq 6000, GRCh37/hg19 alignment, Mutect2-based somatic calling, GATK4 somatic CNV, TruSeq Stranded Total RNA library, STAR alignment, RSEM quantification, FusionCatcher/STAR-Fusion/FusionInspector. Aligned reads, variant calls, segment files and quantifications were obtained by the family and re-analyzed.

**Germline.** Multi-allelic splitting and ClinVar (2026-09-13) annotation with bcftools 1.21; gnomAD v2.1.1 exome frequencies; Ensembl VEP 116 (GRCh37 cache, one transcript per variant, SIFT/PolyPhen, gnomAD MAX_AF). Variants were triaged by ClinVar classification and review status, population frequency, consequence, allele fraction and mapping-quality flags. A 51-gene predisposition panel comprised the ACMG SF v3.2 cancer genes plus sarcoma-associated genes. Germline copy number was screened by per-exon depth (mosdepth) in the normal relative to the sample median.

**Copy number.** snp-pileup on dbSNP common (b151) positions within exon targets; FACETS 0.6.2 [Shen 2016] at cval 150 and 300; whole-genome doubling per Bielski et al. [2018]; HRD components per Telli et al. [2016] re-implemented on exome segments; per-exon tumor:normal depth (mosdepth, normalised to genome-wide exon-target means); rearrangement junctions from split (SA-tag) and discordant reads (samtools).

**Somatic.** VEP 116 annotation of the vendor's somatic VCF; filter as stated in Results; cancer-cell fraction and multiplicity from purity and local total/minor copy number [Dentro 2017]; kataegis screened as ≥3 SNVs within 2 kb.

**Expression.** Tumor RNA re-aligned from FASTQ with STAR 2.7.10b [Dobin 2013] two-pass against GRCh38 with GENCODE v36 using GDC parameters, gene counts converted to TPM with union-exon lengths; chimeric junctions recorded. Cohorts: TARGET-OS and TCGA-SARC STAR-TPM matrices and clinical tables from the UCSC Xena GDC hub [Goldman 2020], primary tumors only. The tumor was quantile-mapped onto the pooled cohort distribution; similarity by Spearman correlation over the 2,000 most variable genes, centroid correlation, and projection onto cohort principal components; per-gene z-scores against each cohort; histone, small-RNA, mitochondrial and pseudoautosomal genes excluded from genome-wide lists.

**Fusions.** FusionCatcher [Nicorici 2014] output scored on unique junction reads, spanning pairs and artifact annotations; STAR chimeric junctions annotated to GENCODE genes as an independent check.

**Immune and HLA.** MCP-counter [Becht 2016] and quanTIseq [Finotello 2019] on the tumor and both cohorts jointly (GDC quantification); OptiType 1.5 [Szolek 2014] on MHC-region reads from normal DNA, tumor DNA and tumor RNA.

**Code and reproducibility.** All scripts, a step-by-step notebook with confidence tags and a corrections log are available at [repository URL]. Compute ran on a single AWS EC2 instance (~4 h total) and a laptop.

---

## Data availability

Processed, de-identified tumor data (somatic variant table, allele-specific copy-number segments, gene-level expression, fusion candidate scores, deconvolution results) will be deposited at [Zenodo DOI]. Raw sequence reads and germline data [will be deposited under controlled access at dbGaP/EGA accession … / are not shared, per the family's consent].

## Ethics statement

[The patient / the patient's legal representative] provided written consent for research use and publication of de-identified genomic data [IRB/ethics statement; some journals accept next-of-kin consent for case reports — confirm with the journal].

## Acknowledgments

[Treating team; Sema4 for data release; TARGET and TCGA contributors; UCSC Xena.]

## Author contributions / Competing interests

[ ]

---

## Figure legends

**Figure 1. Allele-specific copy number.** (A) FACETS genome-wide log-ratio, log-odds-ratio and integer copy-number tracks (cval 150; purity 0.57, ploidy 3.4). *Existing: `results/04_copynumber/facets_cval150.png`.* (B) Chromosome 2 detail: the 18-copy 2p11.2 amplicon with gene annotations (KDM3A, USP39, CYTOR, MAT2A…), the 4→7→18 copy staircase, and the positions of the USP39–CTNNA2 and KDM3A–CYTOR chimeric transcripts. *To make.*

**Figure 2. Biallelic inactivation of TP53 and RB1.** (A) Per-exon tumor:normal log2 depth across TP53 (all exons ≈ −2.1) with the expected value for zero tumor copies at 57% purity. (B) RB1: per-exon depth (exons 1–17 retained, 18+ lost), the p.Thr502Ile variant with allele fraction and inferred multiplicity, and the resolved intron-17 junction (chr13:48,986,411 ↔ 49,237,832, inverted) with split-read support. *To make from `exon_depth_TP53_RB1.tsv` and `rb1_breakpoint.txt`.*

**Figure 3. The MAP4K4 kataegis cluster.** (A) The four SNVs along the gene with strand, trinucleotide context, allele fraction and multiplicity. (B) STAR chimeric reads from MAP4K4 exons to a single downstream acceptor, marking the rearrangement breakpoint relative to the cluster. *To make from `somatic_clean.tsv` and the chimeric junction file.*

**Figure 4. Transcriptome identity.** (A) PCA of TARGET-OS and TCGA-SARC (top 2,000 variable genes) with the ESOS tumor projected; inset: Spearman correlation by cohort group. *Existing: `results/02_expression/gdc_requant/pca_cohorts.png`, `correlation_by_group.png`.* (B) Expression of 2p11.2 amplicon genes vs cohorts. (C) Posterior HOX and osteoblastic genes: tumor vs TARGET-OS vs TCGA-SARC. *To make from `gene_zscores.tsv`, `tumor_extremes_vs_TARGETOS.tsv`.*

**Figure 5. Immune microenvironment.** MCP-counter population scores for the cohorts (boxes) and the tumor (line); quanTIseq fractions as an inset. *Existing: `results/06_hla/immune_mcp.png`.*

**Table 1.** Key genes: integer total/minor copy number, somatic variants, expression percentile vs TARGET-OS. *From `cn_genes.tsv`, `somatic_clean.tsv`, `gene_zscores.tsv`.*

**Supplementary Tables.** S1 FACETS segments; S2 filtered somatic calls with CCF; S3 fusion candidates scored; S4 gene z-scores; S5 deconvolution.

---

## References (to be verified and formatted)

- Becht E et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. *Genome Biol* 2016;17:218.
- Behjati S et al. Recurrent mutation of IGF signalling genes and distinct patterns of genomic rearrangement in osteosarcoma. *Nat Commun* 2017;8:15936.
- Bielski CM et al. Genome doubling shapes the evolution and prognosis of advanced cancers. *Nat Genet* 2018;50:1189–1195.
- Chen X et al. Recurrent somatic structural variations contribute to tumorigenesis in pediatric osteosarcoma. *Cell Rep* 2014;7:104–112.
- Dentro SC, Wedge DC, Van Loo P. Principles of reconstructing the subclonal architecture of cancers. *Cold Spring Harb Perspect Med* 2017;7:a026625.
- Dobin A et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2013;29:15–21.
- Finotello F et al. Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. *Genome Med* 2019;11:34.
- Goldman MJ et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nat Biotechnol* 2020;38:675–678.
- Kovac M et al. Exome sequencing of osteosarcoma reveals mutation signatures reminiscent of BRCA deficiency. *Nat Commun* 2015;6:8940.
- Longhi A et al. Extraskeletal osteosarcoma: a European Musculoskeletal Oncology Society study on 266 patients. *Eur J Cancer* 2017;74:9–16. [verify]
- McLaren W et al. The Ensembl Variant Effect Predictor. *Genome Biol* 2016;17:122.
- Nicorici D et al. FusionCatcher — a tool for finding somatic fusion genes in paired-end RNA-sequencing data. *bioRxiv* 2014. [verify]
- Nik-Zainal S et al. Mutational processes molding the genomes of 21 breast cancers. *Cell* 2012;149:979–993.
- Paludo J et al. Extraskeletal osteosarcoma: outcomes and the role of chemotherapy. *Am J Clin Oncol* 2018;41:832–837. [verify]
- Parrish JK et al. The histone demethylase KDM3A is a microRNA-22-regulated tumor promoter in Ewing sarcoma. *Oncogene* 2015;34:257–262. [verify]
- Perry JA et al. Complementary genomic approaches highlight the PI3K/mTOR pathway as a common vulnerability in osteosarcoma. *PNAS* 2014;111:E5564–E5573.
- Shen R, Seshan VE. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. *Nucleic Acids Res* 2016;44:e131.
- Szolek A et al. OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics* 2014;30:3310–3316.
- Telli ML et al. Homologous recombination deficiency (HRD) score predicts response to platinum-containing neoadjuvant chemotherapy in patients with triple-negative breast cancer. *Clin Cancer Res* 2016;22:3764–3773.
- WHO Classification of Tumours Editorial Board. *Soft Tissue and Bone Tumours*, 5th ed. IARC 2020. [verify chapter]
