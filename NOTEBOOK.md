# Lab notebook — mom_esos

Running log of findings, corrections, open questions and things a human reviewer should check.
Newest entries at the bottom of each section. Dates are when the analysis was run.

**Scope rule:** tumor-level findings are recorded here. Germline variant-level detail is NOT
(README rule: nothing derived from the germline VCF goes into git); it lives only in
`results/01_germline/germline_report.md` on the local disk. This file says only whether the
germline re-check changed anything.

**Confidence tags:** `[solid]` reproduced two ways or expected from first principles;
`[likely]` one method, consistent with everything else; `[weak]` suggestive, needs another method;
`[open]` unresolved.

---

## 1. Status

| Step | Where | State | Report |
|---|---|---|---|
| 00 VEP install | EC2 | done 2026-09-22 (cache: Ensembl 116 GRCh37, 24 GB) | — |
| 01 germline re-check | EC2 | done 2026-09-22, with VEP | `results/01_germline/germline_report.md` |
| 02 expression vs TARGET-OS / TCGA-SARC | Mac | done 2026-09-22; **re-run on GDC-recipe quantification (02b) is the reference version** | `results/02_expression/gdc_requant/expression_report.md` (Sema4-quant run kept alongside) |
| 02b RNA re-quantification (GDC recipe) | EC2 | done 2026-09-22 (163.7 M pairs, 86.8% unique, reverse-stranded) | `results/02_expression/requant/` |
| 03 fusion review | Mac | done 2026-09-22 incl. STAR cross-check; 03b DNA check inconclusive (exome) | `results/03_fusions/fusions_report.md` |
| 04 copy number (FACETS) | EC2 + Mac | done 2026-09-22 | `results/04_copynumber/copynumber_report.md` |
| 05 somatic re-annotation | EC2 + Mac | done 2026-09-22 (VEP + filter policy + CCF); signatures/GENIE comparison pending | `results/05_somatic/somatic_report.md` |
| 06 immune / HLA | EC2 + Mac | done 2026-09-22 (HLA typing; MCP-counter + quanTIseq on the GDC quantification) | `results/06_hla/immune_report.md`, `hla_summary.tsv` (local) |
| 07 deposit + write-up | Mac | not started | |

---

## 2. Findings log

### Case context (family, 2026-09-22; clinical records May–Sept 2020 reviewed 2026-09-22 — source documents are NOT in the repo)
- Primary site: **retroperitoneum**, right upper quadrant, involving pancreatic head, duodenum, IVC (obliterated; bland thrombus below), aorta (near-circumferential), right ureter/kidney (hydronephrosis), right psoas. Unresectable. 12.6 × 10.3 × 7.6 cm at the diagnostic CT (day 0); 17 × 14 cm at day 50 (during RT); 14.2 × 11.9 (day 112) and 13.1 × 8.2 cm (day 143) after RT. No distant metastasis at diagnosis (chest CT, MRI); two indeterminate pulmonary nodules enlarged to 8 and 4 mm at day 143. `[solid]` from radiology reports.
- **Pathology**: initial biopsy day 7 (outside) + CT-guided core day 24 at a tertiary sarcoma center, both high-grade ESOS on expert review. **IHC on the core: SATB2+, SMA multifocal+, desmin scattered+, pan-keratin/AE1-AE3/MDM2/CDK4 negative.** The core biopsy (day 24) is the block Sema4 sequenced → **sequenced tissue is pre-treatment**. `[solid]`
- **Treatment**: primary radiation, 25 fractions over ~6 weeks (days 31–73), interrupted by a 10-day admission at fraction 8 (nausea/vomiting/diarrhea, IVC DVT, small hemoperitoneum). Anticoagulation from day 98. Chemotherapy: reason not documented in the records reviewed. Sema4 report issued day 60. Course after day 143 (Dec 2020 record end) not yet reviewed.
- **Outcome (family, 2026-09-22)**: died of disease ~8.5 months after the diagnostic CT (~day 262), ~6 months after completing radiation. Recorded as an interval only; the date is an identifier and stays out of the repo.
- **Open**: an outside chest/abdomen/pelvis CT ~11 weeks before the diagnostic CT is cited as a comparison in the tertiary-center reports (not available to the community radiologist at diagnosis). Was the mass visible? Request from the tertiary center's outside-imaging archive.
- Implications: (i) DDLPS excluded by IHC as well as genomics; (ii) SMA positivity in *lesional* cells means the MYOCD/ACTA2 transcript signal is likely partly tumor-intrinsic, not only stroma — interpretation revised; (iii) expression data are unaffected by radiation; (iv) an outside CT from ~3 months before diagnosis is referenced in a later report — whether the mass was visible then is an open question for the family.

### Sema4 report, re-read (2026-09-21)
- Report (2020-08-03) calls: homozygous deletions TP53 + RB1; TMB 3.77; MSS; no fusions; 1 VUS (MAP4K4 E589Q); 47-row CN appendix, uninterpreted; "hereditary: none found".
- Report is clinical-scope only. Not done: expression, allele-specific CN / ploidy / LOH, HRD, signatures, immune, cohort comparison. Purity 90% was "per requisition" (pathologist estimate), not computed. `[solid]`
- Sema4 themselves flagged the CN numbers as approximate ("consistent with aneuploidy") and suggested FISH. Not an error on their part.

### 01 Germline (2026-09-22)
- Outcome unchanged from 2020 for cancer predisposition: no P/LP variant in any of 51 panel genes (ACMG cancer genes + sarcoma/OS genes: TP53, RB1, RECQL4, BLM, WRN, CDKN2A, ...). `[solid]` — ClinVar 2026-09-13, gnomAD v2.1.1, VEP 116.
- One incidental carrier-type finding outside the panel, reportable under current (not 2020) practice; for a genetic counselor. Details in the local report only.
- 4 of 6 raw "pathogenic" ClinVar hits were sequencing artifacts (alt fraction 0.17–0.20 at high depth in paralog-rich genes PRSS1, SLC9B1; 5-read GSTT2 call). Sema4's own `RED_FLAGS` (deviantAf, mqTooLow) agreed. Lesson: never read a ClinVar-annotated VCF without allele fraction + mapping-quality flags.
- Germline large-deletion screen (step 04, normal BAM, exon depth vs sample median): no panel gene shows the flat ~0.5x of a heterozygous whole-gene deletion. RB1 0.80x, TP53 unremarkable. POT1 0.64x with patchy exons — reads as capture variability, `[open]`, low prior.

### 04 Copy number (2026-09-22) — FACETS 0.6.2, cval 150 and 300
- **Purity 0.57, ploidy 3.4, whole-genome doubled.** Identical at both cvals. `[solid]` for the qualitative call; purity ±0.05 plausible.
  - Sema4's 90% purity is wrong by a lot. All of Sema4's absolute copy numbers are superseded (diploid baseline + wrong purity).
- **TP53: homozygous deletion of the entire gene** (all 11 exons at log2 ≈ −2.1; predicted for 0 tumor copies at purity 0.57 is ≤ −1.7). 46 kb focal deletion inside a 17p LOH block. `[solid]`
- **RB1: intragenic homozygous deletion starting in intron 17**; exons 1–17 retained at 2 copies with LOH; exons 18–24 at 0 copies `[solid]`; exons 25–27 read ~1 copy on per-exon depth but FACETS calls 0 copies through 49.28 Mb — the three 3' exons are noisy (exon 27 is a 1.9 kb UTR), so "0 copies through exon 27" is `[likely]`.
- **RB1 junction recovered at base-pair resolution (2026-09-22, `scripts/lib/rb1_breakpoint.sh`)**: chr13:48,986,411 (intron 17, hg19) joined to chr13:49,237,832 (~180 kb downstream of RB1's 3' end) in **inverted orientation**; 36 split reads, 94 discordant pairs, MAPQ ≥20. `[solid]` An inverted intrachromosomal junction (fold-back / BFB-type) rather than a simple deletion. Both retained copies carry it (region is LOH at 2 copies), so the event predates WGD. PCR-validatable; output in `results/04_copynumber/rb1_breakpoint.txt`.
- Residual TP53/RB1 RNA (3.7 / 9.9 TPM) is fully explained by the 43% non-tumor cells; earlier suspicion of retained exons withdrawn.
- **KDM3A focal amplification: 18 copies** (Sema4: ×8), 2p11.2, with 96th-percentile expression. `[solid]` amplitude approximate.
- 17p11.2 (AURKB/COPS3) 7 copies + LOH, adjacent to the TP53 deletion — classic OS 17p pattern. 8q (MYC, PRKDC) 6 copies. `[likely]`
- CCNE1, CCND3, RUNX2, VEGFA, CDKN2A: 4 copies = at ploidy, i.e. NOT gained despite Sema4's "gain" labels. `[likely]`
- Single-copy loss + LOH: PTEN, TSC2, CDKN1C. Copy-neutral LOH: ATRX, STK11, CHD5. 32% of autosomal genome is LOH. `[likely]`
- **HRD scar score 75–81** (LOH 19–20, TAI 31–33, LST 25–28; array cutoff 42). `[weak]` — exome segments are coarser than arrays and WGD inflates LST/TAI; but the magnitude is large. BRCA1/2/PALB2 germline clean, so mechanism is open. Compare with Kovac 2015 (BRCAness in ~30% of OS).
- HLA locus (6p21.3) segment: 3 total / 1 minor copy — **no HLA LOH** at segment resolution. `[likely]`; LOHHLA-style allele-specific check not done.

### 06 HLA class I (2026-09-22) — OptiType 1.5 on normal DNA, tumor DNA, tumor RNA
- Heterozygous at HLA-A, -B and -C; **identical 4-digit calls from all three samples**, and all six alleles are expressed in tumor RNA. Combined with the 3/1 FACETS state at 6p21: no HLA loss of any kind. `[solid]` Alleles are in `results/06_hla/hla_summary.tsv` (local only — germline-level information).
- **Immune deconvolution (2026-09-22; MCP-counter + quanTIseq, tumor + both cohorts on the GDC quantification)**: the tumor is *immune-infiltrated for an osteosarcoma* and **myeloid-dominated**. MCP-counter monocytic lineage above every TARGET-OS sample (100th percentile; 97th in TCGA-SARC), neutrophils 94th, NK 91st, T cells 93rd, but CD8 T cells at the OS median (53rd). quanTIseq: macrophages M1 6% + M2 11% (OS medians 1% + 7%), CD4-type T cells 5%, CD8 ~0, dendritic 6%; "Other" 66% is the 2nd percentile of OS, i.e. more non-tumor content than 98% of the cohort — consistent with FACETS purity 0.57. `[likely]` (two methods agree on the myeloid picture; absolute fractions carry the total-RNA caveat).
- Consistent with a large, necrotic retroperitoneal mass (case context). Reading: a macrophage-rich, CD8-poor, PD-L1-low microenvironment — the "myeloid" OS phenotype that does not respond to checkpoint blockade; macrophage-directed or B7-H3-type approaches would be the rational research angle, not PD-1.

### 05 Somatic (2026-09-22) — VEP 116 on Sema4's somatic.vcf, CCF vs FACETS
- **Sema4's PASS set (1,177) is mostly rescued noise.** 1,018 carry `mutectFiltOverride` (Mutect2 filtered them; Sema4's pipeline overrode), 936 `lowAfT`; the flagged calls sit almost entirely below 5% AF and often come as adjacent clusters in one read family (four "ARID1A" calls at 2% AF within 3 codons; four "DICER1"; four "MXI1" frameshifts). `[solid]` These are not mutations. Sema4 reported none of them, so no harm done clinically, but anyone re-using the VCF must filter.
- Filter policy (scripts/05_somatic.R): PASS, ≥8 alt reads, AF ≥0.05 SNV / ≥0.10 indel, normal AF ≤0.02 → **87 calls** (82 SNV, 5 indel). Nonsynonymous TMB from that set 1.8/Mb vs Sema4's 3.77 (their definition unknown). `[likely]`
- **RB1 p.Thr502Ile (c.1505C>T), AF 0.50, 58/117 reads, 0/55 in normal — somatic, clonal, on both retained copies** (CN 2/0 → multiplicity 1.8). COSMIC + HGMD ids; pocket-A domain. **Sema4 did not report it.** So RB1 is inactivated three ways: missense + copy-neutral LOH of exons 1–17 + homozygous deletion of exons 18–24. `[solid]` for the call; pathogenicity of T502I `[likely]` (database-known, pocket domain).
- **MAP4K4 kataegis-like cluster**: 4 clonal missense SNVs within 1.3 kb (E503K, A507T, E517Q, E620Q; Sema4's "E589Q" is one of these under NM_145686 numbering), all G>A/G>C on the + strand = strand-coordinated C>T/C>G, 3 of 4 in TpC context, all at AF 0.47 with 0 normal reads. Multiplicity ~3 of 5 copies → **pre-WGD event**. Mutect2 filtered three of them as clustered events (the known failure mode for kataegis); Sema4 reported one as a VUS. Interpretation: an APOBEC-type localized hypermutation event, not four independent driver mutations. `[likely]` — 4 events is a small cluster for a formal kataegis call.
- No somatic hit in BRCA1/2, PALB2, RAD51C/D, ATM, ATRX, TERT promoter (not captured), TP53 (deleted anyway). HRD mechanism remains `[open]`.
- H3F3A p.Ala115Gly at AF 0.05 (9 reads): not a known hotspot (G34/K27), subclonal at best. Ignore unless it recurs in RNA.
- Substitution spectrum (82 SNVs): C>T 42, C>G 14, C>A 9 — too few for signature fitting beyond "aging-like with some APOBEC"; SigProfiler on this is optional.
- Timing: 24 of 87 kept calls have multiplicity ≥2 in gained regions → pre-WGD; the trunk of this tumor (RB1 hit, MAP4K4 cluster, TP53 loss) predates the doubling. `[likely]`

### Expression, first look (2026-09-21; within-sample percentiles only, no cohort yet)
- Strong osteoblastic program: COL1A1 top gene; SPP1, IBSP, ALPL, RUNX2, SP7 all >93rd percentile. Favors bone-OS-like biology. `[weak]` until step 02.
- T-cell markers near absent (CD8A 0.28 TPM, CD274/PD-L1 0.22), macrophage marker CD163 high, HLA-A/B2M intact. TERT 0.00 TPM (ALT-like telomere maintenance plausible; check ATRX in step 05). CD276/B7-H3 96th percentile.

### 02 Expression vs cohorts (2026-09-22) — TARGET-OS n=88 + TCGA-SARC n=259 (GDC STAR/TPM), tumor quantile-mapped, rank-based
- **The transcriptome is bone-osteosarcoma-like, not soft-tissue-sarcoma-like.** Median Spearman rho (top 2000 variable genes): TARGET-OS 0.62 vs UPS 0.50, MFS 0.48, DDLPS 0.46, LMS 0.39, SS 0.35. 24 of the 25 nearest cohort samples are TARGET-OS. Centroid correlation 0.74 (OS) vs 0.63 (UPS, next). PCA: inside the TARGET-OS cloud, on its UPS/MFS-facing edge (`pca_cohorts.png`). `[solid]` — every method agrees, and a pipeline effect would push the sample away from both cohorts, not toward one.
- Osteoblastic genes are *typical for OS* (z vs OS ≈ 0 for SPP1, IBSP, COL1A1, ALPL, SP7) and extreme vs STS (SP7 z 5.6, IBSP 3.4, SATB2 2.6). RUNX2 is high even for OS (100th percentile) without being amplified (4 copies = ploidy). `[likely]`
- **KDM3A: z +4.8 vs OS, 100th percentile of both cohorts** — the 18-copy focal amplification is expressed. Candidate for the write-up's "novel" line; check literature (KDM3A/JMJD1A, H3K9 demethylase, hypoxia-inducible; Ewing sarcoma link). `[likely]`
- Soft-tissue signal on top of the OS program: MYOCD z 2.9 (3.2 on GDC quant), ACTA2 2.0 (1.2) vs OS (≈ typical for STS; HMGA2 later withdrawn as a pipeline artifact). Clinical IHC (2020) showed **multifocal SMA and scattered desmin positivity in lesional cells**, so the signal is likely partly tumor-intrinsic (myoid differentiation is described in ESOS) as well as retroperitoneal stroma. `[likely]` mixed.
- Immune: PTPRC (CD45) 97th and CD163 97th percentile vs OS, but CD8A/CD3E/NKG7/FOXP3 low (7–30th) and PD-L1 low → macrophage-rich, T-cell-poor. `[likely]`; formal deconvolution pending (06).
- Telomere: TERT 0 TPM; ATRX *high* (z 2.3) → not ATRX-loss ALT; DAXX low (z −2.0 vs OS, −3.1 vs STS). DAXX-loss ALT is a hypothesis worth one line. `[weak]`
- Drug-target genes high vs OS: FGFR1 (98th), PDGFRA (100th), ERBB2 (94th). Research-grade only.
- Odd: E2F1 z −3.3 and CDKN1A (p21) z −2.4 vs OS. p21 low fits TP53 loss; E2F1 low does not fit RB1 loss. Possibly pipeline; `[weak]`.
- Genome-wide extremes on the Sema4 quantification were NOT usable (housekeeping/paralog genes = pipeline differences) → motivated 02b.

### 02 re-run on the GDC-recipe quantification (2026-09-22) — like-for-like with the cohorts; total-RNA vs poly-A remains
- **Sharper, same answer**: median rho TARGET-OS 0.64 vs best STS 0.42 (MFS/UPS); 25/25 nearest neighbours are OS; centroid 0.75 vs 0.52; PCA distance 27 vs 62. `[solid]`
- Key-gene z-scores mostly stable (< 0.5 SD). **Downgraded as pipeline artifacts**: HMGA2 (2.75 → 0.71), VEGFA (0.85 → −0.44), BGLAP, CDK4 (0.13 → −1.19), MDM2 (1.14 → 0.39, i.e. not elevated — no DDLPS-like signal). **Robust**: KDM3A +4.3, RUNX2 +2.3, MYOCD +3.2 (ACTA2 softened to +1.2), PDGFRA +2.0, FGFR1 +2.0, CD163 +2.1, PTPRC +1.6, ATRX +2.1, DAXX −1.8, CDKN1A −2.3, E2F1 −2.4, TP53 −1.4, CDKN1C −2.3, CDK6 +2.2, TOP2A +2.2. ERBB2 dropped to +0.7 (remove from the drug-target line).
- **The whole KDM3A amplicon is transcribed**: six of the top genome-wide "up" genes (TGOLN2, ELMOD3, GGCX, MAT2A, POLR1A, PTCD3; GRCh38 chr2:85.3–86.1 Mb) are amplicon neighbours of KDM3A, plus CYTOR. So the amplicon (~5 Mb, 18 copies) drives coordinate over-expression; KDM3A is the best candidate driver by function, not the only over-expressed gene. 2p11.2 is not a classic OS amplicon (6p12–21, 8q24, 12q13–15, 17p11.2, 19q12 are) → novel. `[likely]`
- **No limb HOX code**: HOXC10 z −4.6, HOXA10 −4.0, HOXA11 −3.1, HOXA9 −2.9, PITX1 −2.5, MSX1 −3.9 vs TARGET-OS (near-zero expression). Skeletal OS carries the limb positional identity of its site of origin; this abdominal tumor does not — same osteoblastic lineage program, different positional identity. `[likely]` (large effects, poly-A mRNAs so library prep is unlikely to explain it). An ESOS-specific observation worth a figure.
- Library-prep signature visible and excluded: histone mRNAs, sn/scaRNA, 7SL, Y_RNA "up" (non-polyA, present in total RNA only); PAR1 genes (CD99, SLC25A6, IL3RA…) "down" because the GENCODE reference carries PAR on X and Y and STAR drops the multimappers. Neither is biology; both are filtered in the script.

### 03 cross-check with the independent STAR alignment (02b) and the DNA check (03b)
- 870 k distant/interchromosomal chimeric reads (total RNA → many rDNA/7SL/collagen artifacts). **No canonical sarcoma fusion** (EWSR1/FUS/SS18/NR4A3/CIC/BCOR/NTRK) in either caller. `[solid]`
- **PPP1CB–ALK (43 reads) is noise**: ALK 0.4 TPM, six different breakpoints, non-canonical junctions — stated explicitly because ALK is targetable.
- CTNNA2–USP39 reproduced by STAR (6 split reads); KDM3A–CYTOR (71 reads) is another intra-amplicon junction. `[likely]` at RNA level.
- **MAP4K4–SLC9A2: 173 STAR split reads sharing one acceptor (GRCh38 chr2:102,453,450) with donors fanning across MAP4K4** → transcription across a genomic breakpoint at the 3' end of MAP4K4, ~10–15 kb from the kataegis cluster. Kataegis adjacent to a rearrangement junction is the expected APOBEC pattern (Nik-Zainal 2012). `[likely]` at RNA level; not seen by FusionCatcher.
- **03b DNA check (exome): inconclusive for all three** (USP39–CTNNA2: 1 tumor pair, 0 normal; PGAP1–DNAH7: 0; MAP4K4 locus: 0). Partner sites are intronic/intergenic, outside capture. Needs WGS or junction PCR. `[open]`

### Fusions, first look (2026-09-21)
- Only FusionInspector-validated call SEC31A--JAK2: 4 junction / 0 spanning reads, non-canonical splice, FFPM 0.02 → artifact. `[solid]` Sema4's "no fusions" stands so far; 605-row raw list not yet reviewed (step 03).

### 03 Fusion review (2026-09-22) — FusionCatcher 176 candidates, Sema4 aggregate 604 rows (= FusionCatcher re-formatted; nothing had a second caller)
- **No oncogenic fusion.** No candidate involves EWSR1/FUS, SS18, NR4A3, CIC, BCOR, NTRK, HEY1–NCOA2 or any other sarcoma-defining partner with more than artifact-level support. Sema4's "none" stands. `[solid]`
- 43 of 176 candidates involve COL1A1 (the tumor's top gene at ~9,000 TPM), with many spanning pairs but ≤8 unique reads: template-switching chimeras of an abundant transcript. Same for NEAT1/MALAT1 lncRNA chimeras.
- **Two "strong" candidates are transcribed DNA breakpoints, not fusions**, and they coincide with FACETS segment boundaries to within a gene: USP39–CTNNA2 (USP39 inside the 18-copy KDM3A amplicon, CTNNA2 at its 4→7-copy flank boundary at 80.10 Mb) and PGAP1–DNAH7 (at the 196.89 / 197.76 Mb 2q segment boundaries). Structural corroboration of the amplicon; PCR-validatable. `[likely]` pending the DNA-level check (03b).

---

### External audit of manuscript v0.1 (2026-09-22) — `manuscript/audit_response_v0.1.md`
- 26-point audit; all its 22 citations verified against PubMed. Accepted essentially all of it; v0.2 written. Analyses added in response: **ASCAT second caller** (04c), **transcriptome robustness** (02c: no-mapping sensitivity, reference-only classifiers with LOO-CV, gene bootstrap, site/age controls, full HOX-cluster table), **RB1 RNA junction analysis** (SJ.out.tab), MAP4K4 transcript reconciliation from the delivered VCF.
- Findings that changed: RB1 p.Thr502Ile downgraded to VUS (gnomAD 1.9e-4, SIFT tolerated; not an inactivating event); "kataegis" → "APOBEC-like localized cluster"; HOX finding is **HOXA-cluster-wide** (+ HOXC9–11), HOXD10/11 elevated, holds vs retroperitoneal STS; CD276 is *below* the OS median (earlier "96th percentile" was within-sample rank); Kovac BRCAness fraction is >80%, not one-third; HRD not "positive", provisional scar burden; HLA "intact antigen presentation" withdrawn (TAP1 below OS median).
- ASCAT vs FACETS: purity 0.67/0.57, ploidy 2.75/3.42, WGD yes/yes, LOH 27%/32%, gene states 9/11 concordant. Reported as ranges. Claim A in REVIEW.md upgraded from "wants a 2nd caller" to "two callers agree qualitatively; third to break the ploidy tie".

## 3. Corrections to earlier statements (kept on purpose)

| Date | Said | Corrected to | Why |
|---|---|---|---|
| 09-21 | TP53 deletion "leaves the last exons" (from Sema4 coords) | Whole gene deleted | Per-exon depth, step 04 |
| 09-21 | "Genome probably doubled; purity may be 55–67%" | WGD confirmed; purity 0.57 | FACETS |
| 09-21 | TP53/RB1 "homozygous deletion unproven, may be partial/subclonal" | Both are clonal homozygous deletions | Log-ratio matches 0 copies at purity 0.57 |
| 09-21 | Residual TP53/RB1 expression "suggests retained exons" | Explained by normal-cell contamination | purity 0.57 |
| 09-22 | RB1 T502I "on both retained copies", one of three inactivating events | Clonal somatic VUS; multiplicity 1.5–1.8 depending on purity; RB1 inactivation rests on the deletion | Audit item 7; ASCAT purity |
| 09-22 | MAP4K4 "kataegis" | "APOBEC-like localized cluster" (4 events do not meet intermutation-distance definitions) | Audit item 11 |
| 09-22 | "Lacks the limb HOX code" / posterior HOX | HOXA-cluster-wide reduction (+HOXC9–11), HOXD10/11 elevated; not a site effect; hypothesis-grade | 02c site-stratified analysis |
| 09-22 | Kovac 2015: BRCAness in "about a third" of OS | >80% (abstract) | Audit item 22 |
| 09-22 | CD276/B7-H3 "96th percentile" | Within-sample rank only; vs TARGET-OS it is at the 19th percentile | 02 GDC z-scores |
| 09-22 | "Antigen presentation intact" | "No HLA class I allele loss detected"; TAP1 below OS median | Audit item 18 |
| 09-22 | "All of Sema4's absolute copy numbers were wrong" | Copy numbers were expressed relative to a diploid reference and flagged approximate by the vendor; allele-specific re-analysis revised them | Audit item 26 |
| 09-22 | MYOCD/ACTA2 signal "stroma, myofibroblastic component less likely" | Partly tumor-intrinsic: clinical IHC shows multifocal SMA in lesional cells | 2020 pathology report |
| 09-22 | "T-cell-poor" microenvironment (from raw CD8A/CD3E TPM) | Immune-infiltrated for OS, myeloid-dominated; CD8 at OS median, CD4-type T cells present | Deconvolution vs cohorts (06) |
| 09-22 | HMGA2 / VEGFA / BGLAP / CDK4 z-scores vs OS | Withdrawn — pipeline artifacts (moved > 1 SD on like-for-like quantification) | 02b |
| 09-22 | Exon-level log2 values in first step-04 run | Values were self-normalised (mosdepth `total_region` covers only the --by BED); fixed to genome-wide exon mean | Bug, see §5 |
| 09-22 | Germline screen "98 exons LOW" | Untargeted UTR/non-coding exons from refGene, not deletions; now judged per gene on captured exons | Method fix |

---

## 4. For expert review (PhD bioinformatician / clinical geneticist)

1. **FACETS solution.** Purity 0.57 / ploidy 3.4 with dipLogR −0.49. FACETS flags "mafR larger than expected if −0.34 is diploid level" — i.e. it considered and rejected a near-diploid solution. Please sanity-check the plots (`results/04_copynumber/facets_cval*.png`); an independent caller (Sequenza or ASCAT) would make this `[solid]`.
2. **HRD score method.** Telli 2016 definitions re-implemented in `scripts/04_copynumber.R` on exome segments; not the scarHRD package. Check the TAI/LST implementations before quoting the number.
3. **RB1 junction** chr13:48,986,411 ↔ 49,237,832, inverted. Please check the read evidence (`rb1_breakpoint.txt`) and the interpretation (fold-back/BFB vs. a two-breakpoint complex event); the exon 25–27 copy state is the ambiguous part. A junction PCR would settle it.
4. **Germline artifact triage** (local report): the calls I dismissed on alt-fraction/mapping grounds (PRSS1 ×2, SLC9B1, GSTT2) — agree?
5. **Germline CNV screen** has no panel of normals; sensitivity for single-exon deletions is poor. POT1 ambiguity.
6. **Somatic filter policy** (step 05, pending): Sema4's PASS set includes calls flagged `readsInN` / `lowAfT`. Decide which flags exclude.
7. **RB1 p.Thr502Ile** — unreported by Sema4 at 50% AF / 117× depth, database-known. Confirm pathogenicity classification (ClinVar/LOVD RB1) and whether it should go on a corrected clinical record.
8. **MAP4K4 cluster = kataegis?** Four strand-coordinated events, 3/4 TpC. Agree with "one APOBEC event" over "four drivers"? Any MAP4K4 sarcoma literature worth citing?
9. **Somatic TMB**: 1.8/Mb by my policy vs Sema4 3.77; which set did Sema4 count? Matters only for the write-up wording.
10. **Expression pipeline mismatch.** The tumor is RSEM/RefSeq, the cohorts are GDC STAR/GENCODE v36. Global similarity is robust (rank-based), single-gene z-scores are indicative, genome-wide differential lists are not. Recommended fix (EC2, ~1 h): re-quantify `ISM563041-2.star.sorted.bam` (or re-align from FASTQ if the BAM's STAR index differs) with the GDC pipeline — STAR 2-pass + GENCODE v36 + HTSeq/STAR counts → TPM — and re-run 02. Also worth checking: was the tumor library poly-A or total RNA (Sema4 lists both kits)? Affects comparability to the poly-A cohorts.
11. **Breakpoint chimeras (03)**: USP39–CTNNA2, PGAP1–DNAH7, KDM3A–CYTOR and MAP4K4–SLC9A2 as transcribed DNA rearrangements — agree? Exome DNA could not confirm (03b). If the block is accessible, a junction PCR on the MAP4K4 3' breakpoint would also anchor the kataegis story.
11b. **HOX/positional identity** (02, GDC run) — needs care. The tumor lacks HOXA10/HOXC10/HOXA9/11/PITX1. The naive reading is "trunk vs limb", but the retroperitoneum is at lumbar level where HOX10 paralogs are normally expressed in axial mesenchyme, so their absence is not explained by site alone and may reflect cell of origin (e.g. a non-skeletal mesenchymal progenitor). Check whether TARGET-OS site metadata lets HOX code be shown to track site within the cohort, and whether TCGA-SARC retroperitoneal tumors (many DDLPS/LMS are retroperitoneal) express HOX10 — that is the right control. Amplicon-wide over-expression: worth a figure.
12. **Sema4 mixed male/female CN files** exist (`segData.female.seg`, `.male.seg`); I used the unsuffixed `segData.seg`. Confirm that is the reported one.

---

## 5. Potential errors, bugs found, and limitations

- **Runner (`run_ec2.sh`)**: SSM ran the remote script under `sh` (pipefail unsupported) — fixed with bash shebang. A `| head -1` inside `set -o pipefail` killed a step silently — fixed. The runner used to continue into the R step after the bash step failed — fixed. It returned exit 0 on remote failure — fixed 09-22.
- **Idle auto-stop alarm** (`esos-idle-stop-*`, 60 min < 3% CPU) stopped the instance mid-download on 09-21 (a download is ~0% CPU). Runner now pauses the alarm's actions for the run and re-arms on exit. Interactive SSH sessions can still trip it.
- **Sema4 VCF quirks**: free-text `## ...` header lines (bcftools rejects), GT written as bare `1` (zygosity is in `INFO/CALL_TYPE`), chrom names `chr*` in VCFs but bare in BAMs. Handled in scripts; anyone re-processing the delivery will hit these.
- **mosdepth `total_region`** is over the `--by` BED only, not genome-wide — bit me once (§3).
- **gnomAD** for the germline step is fetched for panel genes only (remote region read); exome-wide AF comes from VEP's `MAX_AF` (gnomAD exomes r2.1). Two sources, slightly different populations.
- **HRD/WGD from exome**: coarser than array/WGS; treat numbers as approximate.
- **Purity from FACETS vs RNA**: RNA-based deconvolution (step 06) should give a second purity estimate; disagreement would matter.
- **Sema4 somatic VCF `PASS` is not a usable call set** without re-filtering (§2, step 05). Their reportable-variant review clearly used stricter criteria than the VCF FILTER column; the VCF alone would mislead a re-analysis.
- CCF/multiplicity use FACETS cval-150 segments and purity 0.57; if the purity solution moves, every CCF moves with it.
- Expression percentiles so far are within-sample only — not evidence of over-expression until compared with cohorts (step 02).

---

## 6. Run log (EC2)

| Date | Step | Wall time | Notes |
|---|---|---|---|
| 09-21 | 01 | ~1 min ×3 fails, then 1.5 min | sh/pipefail, SIGPIPE, header quirks |
| 09-22 | 00 | ~90 min | first attempt killed by idle alarm at 5 min; second ok |
| 09-22 | 01 (+VEP) | ~4 min | |
| 09-22 | 04 | 13 min first pass (pileup 10 min); 2 short re-runs for normaliser bug | |
| 09-22 | 05 | 1 min | VEP on 1,265 calls |
| 09-22 | 06 | 3 failed starts (conda CLI differs from docs: `optitype run`, razers3 PATH); then 22 min for 3 samples | |
| 09-22 | 04b | 1 min | RB1 breakpoint split-read search |
| 09-22 | 02b | 1 h 45 min | genome download slow (EBI FTP), index 30 min, STAR 2-pass 71 min |
| 09-22 | 03b | 1 min | exome DNA check of RNA chimeras — inconclusive |
| 09-22 | 02c, 04c (Mac) | ~10 min each | audit-response analyses |
| 09-22 | snapshots | — | snap-0956920fd9e1d9979 (root 100 GB), snap-0f1d1d38aa7663a77 (data 1 TB) taken after the batch; instance stopped, volumes kept |

Cost model: instance ~$0.81/h; the 1 TB gp3 volume ~$93/month whether running or not. Plan: finish EC2 batch (05 VEP, 06 HLA), then snapshot/shrink the volume.

---

## 7. Next steps / open questions for the write-up
- ~~Is ESOS-in-this-patient bone-OS-like or soft-tissue-sarcoma-like?~~ **Answered 09-22: bone-OS-like on genome (dual TP53/RB1 loss, WGD, 17p11.2 + 8q amplicons) and transcriptome (step 02).** Remaining nuance: the soft-tissue/myofibroblast signal (MYOCD/ACTA2) — stroma or tumor?
- Mechanism of the high HRD score with clean HR genes germline: somatic hit in BRCA1/2/PALB2/RAD51C? (step 05). ATRX LOH + TERT 0: ALT? (step 05 for an ATRX somatic hit; step 02 for expression).
- Second CN caller for the purity/ploidy solution.
- Deposit plan: processed somatic/CN/expression → open; raw + germline → controlled access, consent needed. PHI in the Sema4 PDF (name, DOB, MRN) never leaves `data/`.
