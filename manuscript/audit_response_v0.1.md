# Response to the v0.1 audit (`case_report_audit_feedback_for_claude.md`, received 2026-09-22)

Manuscript revised to v0.2. Below, each audit point with: our position, what was done, and what remains.
All 22 PMIDs in the audit were checked against PubMed E-utilities on 2026-09-22 and resolved to the stated papers.

| # | Audit point | Position | Done in v0.2 | Open |
|---|---|---|---|---|
| 1 | ESOS treatment literature; add Kondo 2025 | **Accept** | Intro rewritten: chemo benefit supported by retrospective data (Longhi, Paludo, Kondo); regimen class unresolved | — |
| 2 | Prior ESOS genomics under-represented (Jour 2016; Makise 2018) | **Accept** — and Jour 2016 is the paper the Sema4 report itself cited (PMID 27499911) | Both added to Intro and Discussion; novelty claim narrowed to *integrated* WES/WTS + allele-specific CN + breakpoint resolution + cohort classification | — |
| 3 | Validate FACETS with a second/third caller | **Accept** | **ASCAT run on the same allele counts** (`scripts/04c_ascat.R`, `results/04_copynumber/ascat_report.md`): purity 0.67 / ploidy 2.75 vs FACETS 0.57 / 3.42; single likelihood optimum, no near-diploid solution; WGD yes in both; LOH 27% vs 32%; integer states agree at 9/11 loci (TP53 46-kb deletion below ASCAT segment resolution; KDM3A 16 vs 18). Manuscript now reports ranges. | Third caller (Sequenza or PureCN, needs EC2) to break the ploidy tie; concordance table template in `ascat_vs_facets.tsv` |
| 4 | "Sema4's 90% purity was wrong" too strong; "assumed diploid baseline" undemonstrated | **Accept the wording; keep the diploid-reference statement, now evidenced** | "Requisition tumor-cellularity estimate was 90%; allele-specific analysis inferred 0.57–0.67." The diploid reference is demonstrable from the delivered segment file (`callType == normal` rows carry `numCopies ≈ 2.0`) and is stated that way | — |
| 5 | KDM3A: prefer OS-specific literature; "novel" → "apparently non-recurrent" | **Accept** | Wang 2022 and Li 2022 added; Parrish 2015 kept as secondary; wording changed; clinical report credited for detecting the amplicon | Systematic check of 2p11.2 gain frequency in TARGET-OS CN data (Xena GDC masked-CNV) — not done |
| 6 | TP53/RB1 structural result strong | Noted | Kept as major result | — |
| 7 | RB1 p.Thr502Ile over-interpreted | **Accept** | Reclassified as "clonal somatic variant of uncertain functional significance" (gnomAD 1.9×10⁻⁴, SIFT tolerated, PolyPhen possibly damaging, conflicting ClinVar); no longer counted as an inactivating event; multiplicity given as 1.5–1.8 across the two purity estimates ("consistent with, not proven") | ClinVar-Miner/LOVD classification to be cited explicitly |
| 8 | Timing of deletions vs WGD needs event-specific evidence | **Partly accept** | The event-specific evidence exists and is now stated: the RB1 5' segment is at 2 copies with LOH while the 3' segment is at 0, so both retained 13q copies carry the same junction → deletion before arm duplication (most parsimonious). Same logic for TP53 within its 17p LOH block, labelled "less certain" | — |
| 9 | RB1 RNA attribution: exons 1–17 are retained and may be transcribed | **Accept — good catch** | RNA splice-junction analysis added (`SJ.out.tab`): junction reads across exons 18–27 are 52% of those across exons 1–17, vs 43% expected if 3' transcripts come only from non-tumor cells; reads splice from RB1 exons into the intron-17 rearrangement region. Text now says 5' RB1 is transcribed in tumor cells | Per-exon RNA coverage plot (Supp Fig S2) needs the RNA BAM, which was deleted to save space — re-run 02b with BAM retained, or use the Sema4 hg19 STAR BAM on EC2 |
| 10 | TMB 1.8 vs 3.77 not comparable | **Accept** | Our definition stated explicitly (non-synonymous filtered calls / 35 Mb); comparison framed as "not directly comparable"; no claim of which is correct | Reproducing the vendor definition is not possible (unpublished) |
| 11 | "Kataegis" stronger than four mutations support | **Accept** | "APOBEC-like localized mutation cluster"; text acknowledges the intermutation-distance definitions are not met; in-cis phasing not possible with exome read length (events 12–1,300 bp apart; the two closest, 102,476,129 and 102,476,141, could be phased on single reads — see open) | Phasing of the two 12-bp-apart SNVs from read pairs (EC2, trivial); exact breakpoint needs WGS |
| 12 | MAP4K4 transcript reconciliation must be auditable | **Accept** | Genomic positions given in text; Sema4's vkey `_266uu066uu0012` (c.1765G>C, NM_145686.3, p.Glu589Gln) = chr2:102,477,440 G>C = our p.Glu620Gln (ENST00000314363), verified from the delivered VCF/annotation files; Supplementary Table S3 specified | Build S3 as a file |
| 13 | "No gene fusion was present" too absolute; ETV6::NTRK3 ESOS case | **Accept** | "No convincing oncogenic or sarcoma-defining fusion was detected"; Skok 2025 cited in Intro | — |
| 14–15 | Harden transcriptome classification | **Accept** | `scripts/02c_robustness.R` → `robustness_report.md`: (i) raw vs quantile-mapped: Spearman identical, PCA distances near-identical; (ii) reference-only nearest-centroid and kNN, LOO-CV 99.4% / 99.1% for OS vs STS, ESOS → OS by both (15/15 votes); (iii) 500 gene-resampling draws: OS top in 100%, margin 0.15 (0.12–0.18); (iv) cohort age/site distributions reported; the confound works against the OS call | Independent (adult) OS expression cohort; axial/pelvic OS is n=5 in TARGET — too small |
| 16 | HOX interpretation too strong | **Accept — and the analysis changed the finding** | Full four-cluster table + site-stratified comparison: the reduction is **HOXA-cluster-wide** (every HOXA gene below OS mean; A9/A10/A11 at −2.7 to −3.4) plus HOXC9–11, with HOXB unremarkable and HOXD10/11 *elevated*; it holds against retroperitoneal TCGA-SARC (A9/A10/A11 −1.4 to −2.1). Rewritten as "reduced HOXA-cluster expression"; cell-of-origin and epigenetic silencing presented as hypotheses; age matching and HOXA methylation listed as needed | Age-adjusted comparison (no adult OS cohort in hand); HOXA locus methylation |
| 17 | Deconvolution ≠ cell counts | **Accept** | "quanTIseq estimated a macrophage fraction of approximately 17%"; IHC panel (CD68, CD163, CD8, B7-H3) listed in Limitations | IHC if block accessible |
| 18 | "Antigen presentation intact" too strong | **Accept** | "No HLA class I allele loss was detected; all six alleles had RNA support"; TAP1 below OS median now stated | — |
| 19 | Checkpoint discussion | **Accept** | SARC028 / PEMBROSARC cited; no claim that this tumor's profile predicts non-response | — |
| 20 | B7-H3 | **Accept — and corrected a fact** | CD276 is expressed (57 TPM) but *below the osteosarcoma median* (19th percentile vs OS); the earlier "96th percentile" was within-sample rank. Discussion says so; trials kept generic | — |
| 21 | HRD score interpretation | **Accept** | Not labelled HRD-positive; "high provisional genomic-scar burden"; Sztupinszki 2018 and Nacer 2026 cited; validation path (scarHRD on Sequenza output; HRDetect/CHORD with WGS; position within an OS cohort) listed | scarHRD requires Sequenza (EC2); OS-cohort comparison needs allele-specific segments for TARGET-OS |
| 22 | Kovac "one-third" wrong | **Accept — our error** | Corrected to ">80%" (verified from the abstract) | — |
| 23 | PARP as research implication only | **Accept** | Done; Engert 2017 cited | — |
| 24 | Germline wording; ACMG v3.3 | **Accept** | "Research reanalysis detected no P/LP coding/splice variant in the prespecified genes"; exclusions listed; panel labelled v3.2 with a note to check v3.3 | Update panel to v3.3 and re-run step 01 (minutes) |
| 25 | Generalization from n=1 | **Accept** | Title, abstract and discussion reframed as "striking convergence with conventional osteosarcoma" for this case; heterogeneity literature cited | — |
| 26 | Tone toward the clinical report | **Accept** | Neutral wording throughout; clinical report credited for TP53/RB1, the 2p11.2/KDM3A amplification and its aneuploidy caveat; the VCF re-filtering point kept as a factual data note in Methods/Discussion | — |

## Where we push back (mildly)

- **Spearman and quantile mapping (14/15):** the mapping is a monotone transform of the tumor vector, so Spearman is unaffected by construction; the sensitivity analysis confirms identical values. The PCA result was the one that could have moved, and it did not (OS distance 26.7 raw vs 27.1 mapped).
- **Deletion timing (8):** we agree "final copy state alone" is weak in general, but here the 5'/3' RB1 copy-state contrast is event-specific evidence that both copies carry the junction; we state it as parsimony, not proof.
- **Diploid reference (4):** the vendor's convention is visible in their own segment file, so we keep the statement, worded descriptively.

## Priority list after this revision (in the order we would do them)

1. Third CN caller on EC2 (Sequenza → also feeds scarHRD) — ~1 h instance time.
2. Re-run step 01 with the ACMG SF v3.3 gene set.
3. Phase the two closest MAP4K4 SNVs from read pairs (EC2, minutes); build Supplementary Table S3.
4. Per-exon RB1 RNA coverage figure (needs the RNA BAM retained in a 02b re-run).
5. Pathology/IHC and junction PCR on the archived block, if accessible and consented.
6. Independent adult osteosarcoma expression cohort for age-matched comparison.
