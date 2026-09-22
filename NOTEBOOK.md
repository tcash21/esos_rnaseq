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
| 02 expression vs TARGET-OS / TCGA-SARC | Mac | not started | |
| 03 fusion review | Mac | not started (one FusionInspector call already judged artifact, see §2) | |
| 04 copy number (FACETS) | EC2 + Mac | done 2026-09-22 | `results/04_copynumber/copynumber_report.md` |
| 05 somatic re-annotation | EC2 + Mac | done 2026-09-22 (VEP + filter policy + CCF); signatures/GENIE comparison pending | `results/05_somatic/somatic_report.md` |
| 06 immune / HLA | EC2 + Mac | HLA typing done 2026-09-22; immune deconvolution not started | `results/06_hla/hla_summary.tsv` |
| 07 deposit + write-up | Mac | not started | |

---

## 2. Findings log

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
- Immune deconvolution (quanTIseq / MCP-counter on genes.results) is the Mac part, not started.

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

### Fusions, first look (2026-09-21)
- Only FusionInspector-validated call SEC31A--JAK2: 4 junction / 0 spanning reads, non-canonical splice, FFPM 0.02 → artifact. `[solid]` Sema4's "no fusions" stands so far; 605-row raw list not yet reviewed (step 03).

---

## 3. Corrections to earlier statements (kept on purpose)

| Date | Said | Corrected to | Why |
|---|---|---|---|
| 09-21 | TP53 deletion "leaves the last exons" (from Sema4 coords) | Whole gene deleted | Per-exon depth, step 04 |
| 09-21 | "Genome probably doubled; purity may be 55–67%" | WGD confirmed; purity 0.57 | FACETS |
| 09-21 | TP53/RB1 "homozygous deletion unproven, may be partial/subclonal" | Both are clonal homozygous deletions | Log-ratio matches 0 copies at purity 0.57 |
| 09-21 | Residual TP53/RB1 expression "suggests retained exons" | Explained by normal-cell contamination | purity 0.57 |
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
10. **Sema4 mixed male/female CN files** exist (`segData.female.seg`, `.male.seg`); I used the unsuffixed `segData.seg`. Confirm that is the reported one.

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

Cost model: instance ~$0.81/h; the 1 TB gp3 volume ~$93/month whether running or not. Plan: finish EC2 batch (05 VEP, 06 HLA), then snapshot/shrink the volume.

---

## 7. Next steps / open questions for the write-up
- Is ESOS-in-this-patient bone-OS-like or soft-tissue-sarcoma-like? Current evidence (dual TP53/RB1 loss, WGD, 17p11.2 + 8q + KDM3A amplification, osteoblastic expression) points to bone-OS. Step 02 cohort comparison is the test.
- Mechanism of the high HRD score with clean HR genes germline: somatic hit in BRCA1/2/PALB2/RAD51C? (step 05). ATRX LOH + TERT 0: ALT? (step 05 for an ATRX somatic hit; step 02 for expression).
- Second CN caller for the purity/ploidy solution.
- Deposit plan: processed somatic/CN/expression → open; raw + germline → controlled access, consent needed. PHI in the Sema4 PDF (name, DOB, MRN) never leaves `data/`.
