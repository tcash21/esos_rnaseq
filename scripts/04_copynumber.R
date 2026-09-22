#!/usr/bin/env Rscript
# =============================================================================
# 04_copynumber.R — interpret the FACETS fit (runs on the Mac; inputs are small tables).
#   Rscript scripts/04_copynumber.R            (after: bash scripts/run_ec2.sh 04)
#
# Reads results/04_copynumber/facets_cval*_{fit,segments}.tsv, exon_depth_TP53_RB1.tsv,
# normal_exon_depth_panel.tsv, plus Sema4's GATK .seg from data/sema4.
# Writes results/04_copynumber/
#   cn_genes.tsv              integer total/minor CN for OS-relevant genes, per cval, vs Sema4's call
#   cn_arms.tsv               arm-level total CN (for TARGET-OS comparison in the write-up)
#   hrd_scores.tsv            LOH / TAI / LST counts per cval (Telli 2016 definitions, exome-approximate)
#   copynumber_report.md      purity, ploidy, WGD, LOH fraction, TP53/RB1 exon map, germline CNV screen
# =============================================================================
suppressPackageStartupMessages(library(data.table))
proj <- Sys.getenv("PROJ", getwd())
out  <- file.path(proj, "results/04_copynumber")
rd   <- function(f, ...) fread(file.path(out, f), ...)

fits <- rbindlist(lapply(list.files(out, "^facets_cval[0-9]+_fit\\.tsv$", full.names = TRUE), fread))
setorder(fits, cval)
segs <- rbindlist(lapply(fits$cval, function(cv) rd(sprintf("facets_cval%d_segments.tsv", cv))[, cval := cv]))
segs[, `:=`(len = end - start, chrom = as.character(chrom))]
segs[, lcn.em := as.numeric(lcn.em)]                     # NA where FACETS could not resolve the minor allele

# ---------- genome-level calls ----------
# WGD (Bielski 2018): >50% of the autosomal genome has major CN >= 2
gl <- segs[!chrom %in% c("X", "Y") & !is.na(tcn.em), .(
  ploidy_wtd   = sum(tcn.em * len) / sum(len),
  frac_majorCN_ge2 = sum(len[ifelse(is.na(lcn.em), tcn.em, tcn.em - lcn.em) >= 2]) / sum(len),
  frac_LOH     = sum(len[!is.na(lcn.em) & lcn.em == 0]) / sum(len),
  frac_gain    = sum(len[tcn.em > round(sum(tcn.em * len) / sum(len))]) / sum(len),
  frac_loss    = sum(len[tcn.em < round(sum(tcn.em * len) / sum(len))]) / sum(len),
  frac_homdel  = sum(len[tcn.em == 0]) / sum(len)), by = cval]
gl[, WGD := frac_majorCN_ge2 > 0.5]
gl <- merge(fits[, .(cval, purity, ploidy, dipLogR, flags)], gl, by = "cval")

# ---------- HRD scores (Telli 2016; exome-approximate — segment boundaries are coarser than SNP arrays) ----------
hrd <- segs[!chrom %in% c("X", "Y") & !is.na(lcn.em), {
  chr_len <- .SD[, .(L = max(end) - min(start)), by = chrom]
  # LOH: lcn==0 segments > 15 Mb that are not the whole chromosome
  loh <- .SD[lcn.em == 0 & len > 15e6][chr_len, on = "chrom"][len < 0.9 * L, .N]
  # TAI: allelic imbalance (lcn != tcn - lcn) extending to a telomere, not crossing the centromere (approx: segment touches chrom end)
  tai <- .SD[, { s <- .SD[order(start)]; imb <- s$tcn.em - 2 * s$lcn.em != 0
                 sum(imb[c(1, .N)]) }, by = chrom][, sum(V1)]
  # LST: breakpoints between adjacent segments both >= 10 Mb (after dropping < 3 Mb segments)
  lst <- .SD[len >= 3e6][order(chrom, start), { s <- .SD; if (.N < 2) 0L else
                 sum(s$len[-1] >= 10e6 & s$len[-.N] >= 10e6 & (s$tcn.em[-1] != s$tcn.em[-.N] | s$lcn.em[-1] != s$lcn.em[-.N])) }, by = chrom][, sum(V1)]
  .(LOH = loh, TAI = tai, LST = lst, HRD_sum = loh + tai + lst)
}, by = cval]

# ---------- gene-level ----------
genes <- fread(text = "gene\tchrom\tstart\tend
TP53\t17\t7571720\t7590868
RB1\t13\t48877887\t49056122
CDKN2A\t9\t21967751\t21995300
MDM2\t12\t69201952\t69239214
CDK4\t12\t58141510\t58146230
MYC\t8\t128748315\t128753680
RUNX2\t6\t45328330\t45551082
VEGFA\t6\t43737946\t43754224
CCND3\t6\t41902671\t42018262
CCNE1\t19\t30302805\t30315215
KDM3A\t2\t86668435\t86719405
PTEN\t10\t89623195\t89728532
ATRX\tX\t76760356\t77041719
DLG2\t11\t83166055\t85338314
AKT3\t1\t243651535\t244006886
CHD5\t1\t6161846\t6240194
AURKB\t17\t8108042\t8113990
COPS3\t17\t17151024\t17169365
PDGFRA\t4\t55095264\t55164412
IGF1R\t15\t99192200\t99507759
STK11\t19\t1205798\t1228434
TSC2\t16\t2097990\t2139492
CDKN1C\t11\t2904443\t2907111
NF1\t17\t29421945\t29704695
PRKDC\t8\t48685969\t48872743")
ov <- segs[genes, on = .(chrom, start <= end, end >= start), nomatch = NULL,
           .(gene, cval, tcn = tcn.em, lcn = lcn.em, cf = cf.em, cnlr = cnlr.median, seg_len = len)]
# a gene spanning a breakpoint gets the segment with the LOWER tcn (deletion is what we care about)
cn_genes <- ov[order(gene, cval, tcn)][, .SD[1], by = .(gene, cval)]
cn_genes[, state := fcase(tcn == 0, "HOMDEL", !is.na(lcn) & lcn == 0 & tcn == 1, "HET_LOSS+LOH",
                          !is.na(lcn) & lcn == 0, "LOH", tcn >= 2 * round(gl[cval == .BY$cval, ploidy_wtd]) + 2, "AMP",
                          default = ""), by = cval]

# Sema4's GATK numCopies for the same genes
s4f <- list.files(file.path(proj, "data/sema4"), "cnvSomatic\\.segData\\.seg$", recursive = TRUE, full.names = TRUE)[1]
if (!is.na(s4f)) {
  s4 <- fread(s4f)[sample == "ISM556046-2", .(chrom = as.character(segChrom), start = segL, end = segR, sema4_copies = numCopies)]
  s4g <- s4[genes, on = .(chrom, start <= end, end >= start), nomatch = NULL, .(gene, sema4_copies)][order(gene, sema4_copies)][, .SD[1], by = gene]
  cn_genes <- merge(cn_genes, s4g, by = "gene", all.x = TRUE)
}
fwrite(cn_genes[order(cval, gene)], file.path(out, "cn_genes.tsv"), sep = "\t")

# ---------- arm-level ----------
cen <- fread(text = "chrom\tcen
1\t125000000\n2\t93300000\n3\t91000000\n4\t50400000\n5\t48400000\n6\t61000000\n7\t59900000\n8\t45600000\n9\t49000000\n10\t40200000
11\t53700000\n12\t35800000\n13\t17900000\n14\t17600000\n15\t19000000\n16\t36600000\n17\t24000000\n18\t17200000\n19\t26500000\n20\t27500000
21\t13200000\n22\t14700000\nX\t60600000")
arms <- segs[cen, on = "chrom", nomatch = NULL][, arm := ifelse(end <= cen, "p", ifelse(start >= cen, "q", "pq"))][arm != "pq",
        .(tcn_wtd = sum(tcn.em * len) / sum(len), frac_LOH = sum(len[!is.na(lcn.em) & lcn.em == 0]) / sum(len), n_seg = .N),
        by = .(cval, chrom, arm)]
fwrite(arms[order(cval, as.integer(chrom), arm)], file.path(out, "cn_arms.tsv"), sep = "\t")
fwrite(hrd, file.path(out, "hrd_scores.tsv"), sep = "\t")

# ---------- TP53 / RB1 exon map & germline screen ----------
ex <- if (file.exists(file.path(out, "exon_depth_TP53_RB1.tsv"))) rd("exon_depth_TP53_RB1.tsv") else NULL
gs <- if (file.exists(file.path(out, "normal_exon_depth_panel.tsv"))) rd("normal_exon_depth_panel.tsv") else NULL

# ---------- report ----------
f2 <- function(x) sprintf("%.2f", x); pc <- function(x) sprintf("%.0f%%", 100 * x)
exon_line <- function(g) if (is.null(ex)) "_not run_" else ex[grepl(paste0("^", g, "_"), exon)][order(as.integer(sub(".*exon", "", exon))),
  paste0(sub(".*_", "", exon), ":", ifelse(is.na(log2_ratio), "NA", f2(log2_ratio)))] |> paste(collapse = "  ")
rep <- c(
  "# Copy number — FACETS allele-specific fit vs Sema4 GATK calls",
  paste0("_generated ", Sys.Date(), "; hg19; research-grade_"), "",
  "## Purity / ploidy (one row per segmentation strictness; agreement between rows = stable solution)",
  "| cval | purity | ploidy | dipLogR | WGD | genome LOH | gain | loss | homdel | flags |", "|---|---|---|---|---|---|---|---|---|---|",
  gl[, sprintf("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |", cval, f2(purity), f2(ploidy), f2(dipLogR), ifelse(WGD, "**yes**", "no"),
               pc(frac_LOH), pc(frac_gain), pc(frac_loss), pc(frac_homdel), ifelse(flags == "", "-", flags))], "",
  "Sema4 assumed purity 90% (per requisition) and a diploid baseline. WGD = >50% of autosomes with major allele CN >= 2 (Bielski 2018).", "",
  "## HRD scar scores (Telli 2016 definitions, exome-approximate; array-based cutoff HRD_sum >= 42)",
  "| cval | LOH | TAI | LST | sum |", "|---|---|---|---|---|",
  hrd[, sprintf("| %d | %d | %d | %d | **%d** |", cval, LOH, TAI, LST, HRD_sum)], "",
  "## Key genes (integer total / minor copy number; lcn NA = unresolved)",
  "| gene | cval | tcn | lcn | state | cell frac | Sema4 copies |", "|---|---|---|---|---|---|---|",
  cn_genes[order(gene, cval), sprintf("| %s | %d | %d | %s | %s | %s | %s |", gene, cval, tcn, ifelse(is.na(lcn), "NA", lcn), state,
                                      ifelse(is.na(cf), "NA", f2(cf)), if ("sema4_copies" %in% names(cn_genes)) f2(sema4_copies) else "-")], "",
  "## TP53 / RB1 per-exon tumor:normal depth (log2; ~0 = same as normal, -1 = half, < -2 = near-absent)",
  paste0("- **TP53** ", exon_line("TP53")), paste0("- **RB1** ", exon_line("RB1")), "",
  "## Germline CNV screen (normal BAM, predisposition panel exons vs sample median; no panel-of-normals so coarse)",
  if (is.null(gs)) "_not run_" else {
    # refGene lists UTR/non-coding exons the capture never targeted (depth ~0), so judge per gene on captured exons only.
    # A heterozygous whole-gene deletion = every captured exon near 0.5x; a single-exon deletion = one exon near 0.5x with good depth.
    pg <- gs[normal_depth >= 20, .(n_exons = .N, median_rel = median(rel_to_median), n_below_0.7 = sum(rel_to_median < 0.7)), by = gene][order(median_rel)]
    sus <- pg[median_rel < 0.65]
    c(sprintf("- %d panel genes assessed on captured exons; lowest per-gene median relative depth: %s",
              nrow(pg), pg[1:5, paste0(gene, " ", sprintf("%.2f", median_rel), collapse = ", ")]),
      if (nrow(sus)) sus[, sprintf("- **%s** median %.2fx over %d exons (%d below 0.7x) — possible germline deletion or poor capture; check het SNP allele fractions across the gene", gene, median_rel, n_exons, n_below_0.7)]
      else "- No gene has its captured exons consistently near 0.5x: no evidence of a heterozygous germline deletion of TP53, RB1 or any other panel gene (sensitivity for single-exon events is limited without a panel of normals)")
  }, "",
  "## Reading guide",
  "- If purity comes out well under 0.9, Sema4's 'copies = 0.6' for TP53/RB1 is consistent with a true homozygous deletion diluted by normal cells; if purity is ~0.9, the deletion is subclonal or partial — see the exon map.",
  "- Arm-level table (cn_arms.tsv) is what to compare against TARGET-OS / Sema4 appendix in the write-up.")
writeLines(rep, file.path(out, "copynumber_report.md")); cat(rep, sep = "\n")
