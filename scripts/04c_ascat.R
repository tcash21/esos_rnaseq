#!/usr/bin/env Rscript
# =============================================================================
# 04c_ascat.R — independent purity/ploidy/allele-specific CN with ASCAT, from the same snp-pileup as FACETS (Mac).
#   PROJ=$PWD Rscript scripts/04c_ascat.R
#
# Audit v0.1 priority 1: reproduce the FACETS solution with a second caller. ASCAT (Van Loo 2010) models LogR and
# BAF jointly with its own purity/ploidy grid search — a different statistical route to the same quantities.
# Input : results/04_copynumber/facets_pileup.csv.gz  (File1 = normal, File2 = tumor; R/A ref/alt counts)
# Output: results/04_copynumber/ascat/  (ASCAT plots, segments) + ascat_vs_facets.tsv, ascat_report.md
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ASCAT) })
proj <- Sys.getenv("PROJ", getwd()); cn <- file.path(proj, "results/04_copynumber"); out <- file.path(cn, "ascat")
dir.create(out, showWarnings = FALSE, recursive = TRUE); setwd(out)

p <- fread(cmd = paste("gzcat", shQuote(file.path(cn, "facets_pileup.csv.gz"))))
p <- p[Ref != "." & Chromosome %in% c(1:22, "X")]
p[, `:=`(nd = File1R + File1A, td = File2R + File2A)]
p <- p[nd >= 25 & td >= 25]                                   # same depth floor as the FACETS run
p[, `:=`(baf_n = File1A / nd, baf_t = File2A / td)]
# LogR: tumour/normal depth ratio, median-centred (ASCAT re-centres via ploidy anyway); BAF as alt fraction.
p[, logr := log2((td / sum(td)) / (nd / sum(nd)))]
p[, logr := logr - median(logr)]
# ASCAT expects germline genotypes to pick heterozygous SNPs: call hets in the normal at 0.3–0.7 with >= 25 reads.
p[, het := baf_n > 0.3 & baf_n < 0.7]
p[, id := paste0("snp", .I)]
cat(sprintf("positions after depth filter: %d; germline hets: %d\n", nrow(p), sum(p$het)))

w <- function(v) { m <- matrix(v, ncol = 1, dimnames = list(p$id, "ESOS")); m }
Tumor_LogR <- w(p$logr); Tumor_BAF <- w(p$baf_t); Germline_LogR <- w(0 * p$logr); Germline_BAF <- w(p$baf_n)
snppos <- data.frame(chrs = p$Chromosome, pos = p$Position, row.names = p$id)
for (nm in c("Tumor_LogR", "Tumor_BAF", "Germline_LogR", "Germline_BAF")) {
  d <- cbind(snppos, get(nm)); write.table(d, paste0(nm, ".txt"), sep = "\t", quote = FALSE, col.names = NA) }
asc <- ascat.loadData("Tumor_LogR.txt", "Tumor_BAF.txt", "Germline_LogR.txt", "Germline_BAF.txt", chrs = c(1:22, "X"), gender = "XX", genomeVersion = "hg19")
ascat.plotRawData(asc, img.prefix = "raw_")
gg <- list(germlinegenotypes = matrix(!p$het, ncol = 1, dimnames = list(p$id, "ESOS")), failedarrays = NULL)
asc <- ascat.aspcf(asc, ascat.gg = gg, penalty = 70)         # penalty ~ exome default range
ascat.plotSegmentedData(asc, img.prefix = "seg_")
res <- ascat.runAscat(asc, gamma = 1, img.prefix = "ascat_", write_segments = TRUE)   # gamma = 1 for sequencing data
seg <- as.data.table(res$segments)
setnames(seg, c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor"))
seg[, `:=`(len = endpos - startpos, tcn = nMajor + nMinor)]
fwrite(seg, file.path(cn, "ascat_segments.tsv"), sep = "\t")

# ---------- comparison with FACETS ----------
fa <- fread(file.path(cn, "facets_cval150_fit.tsv")); fs <- fread(file.path(cn, "facets_cval150_segments.tsv"))[, chrom := as.character(chrom)]
auto <- seg[chr != "X"]
wgd_ascat <- auto[, sum(len[nMajor >= 2]) / sum(len)] > 0.5
loh_ascat <- auto[, sum(len[nMinor == 0]) / sum(len)]
fs_auto <- fs[chrom != "X"]; fs_auto[, len := end - start]
loh_facets <- fs_auto[, sum(len[!is.na(lcn.em) & lcn.em == 0]) / sum(len)]
genes <- fread(text = "gene\tchrom\tstart\tend
TP53\t17\t7571720\t7590868
RB1_5prime\t13\t48877887\t48986400
RB1_3prime\t13\t48986500\t49056122
KDM3A\t2\t86668435\t86719405
AURKB_17p11\t17\t8108042\t8113990
MYC\t8\t128748315\t128753680
PTEN\t10\t89623195\t89728532
ATRX\tX\t76760356\t77041719
CDKN2A\t9\t21967751\t21995300
RUNX2\t6\t45328330\t45551082
MAP4K4\t2\t102313000\t102511000")
pick <- function(d, sc, ss, se, tc, lc) { r <- d[get(sc) == genes$chrom[i] & get(ss) <= genes$end[i] & get(se) >= genes$start[i]]
  if (nrow(r) == 0) return(c(NA, NA)); r <- r[order(get(tc))][1]; c(r[[tc]], r[[lc]]) }
cmp <- rbindlist(lapply(seq_len(nrow(genes)), function(j) { i <<- j
  f <- pick(fs, "chrom", "start", "end", "tcn.em", "lcn.em"); a <- pick(seg[, chr := as.character(chr)], "chr", "startpos", "endpos", "tcn", "nMinor")
  data.table(gene = genes$gene[j], facets_tcn = f[1], facets_lcn = f[2], ascat_tcn = a[1], ascat_lcn = a[2]) }))
summary_dt <- data.table(metric = c("purity", "ploidy", "WGD", "genome LOH fraction", "n segments"),
                         FACETS = c(round(fa$purity, 3), round(fa$ploidy, 2), "yes", round(loh_facets, 3), nrow(fs)),
                         ASCAT = c(round(res$aberrantcellfraction, 3), round(res$ploidy, 2), ifelse(wgd_ascat, "yes", "no"), round(loh_ascat, 3), nrow(seg)))
fwrite(rbind(summary_dt, data.table(metric = paste0(cmp$gene, " tcn/lcn"), FACETS = paste0(cmp$facets_tcn, "/", cmp$facets_lcn), ASCAT = paste0(cmp$ascat_tcn, "/", cmp$ascat_lcn))),
       file.path(cn, "ascat_vs_facets.tsv"), sep = "\t")
rep <- c("# ASCAT vs FACETS — independent purity / ploidy / allele-specific copy number", paste0("_generated ", Sys.Date(), "; same snp-pileup input; ASCAT ", as.character(packageVersion("ASCAT")), ", gamma 1, ASPCF penalty 70_"), "",
  "| metric | FACETS (cval 150) | ASCAT |", "|---|---|---|", summary_dt[, sprintf("| %s | %s | %s |", metric, FACETS, ASCAT)], "",
  "| locus | FACETS total/minor | ASCAT total/minor |", "|---|---|---|", cmp[, sprintf("| %s | %s/%s | %s/%s |", gene, facets_tcn, facets_lcn, ascat_tcn, ascat_lcn)], "",
  sprintf("ASCAT goodness of fit: %.1f%%. Plots in results/04_copynumber/ascat/ (ascat_ESOS.png shows the purity/ploidy sunrise plot).", res$goodnessOfFit),
  "",
  "Reading (2026-09-22 run): both callers place the tumour in the same purity/ploidy basin (mid purity, hyperdiploid, WGD) and neither offers a near-diploid / ~90%-purity solution (see the sunrise plot). Exact values differ along the purity–ploidy trade-off (ASCAT assigns more of the genome to 3 copies where FACETS assigns 4). Integer states agree at 9/11 loci; the two exceptions are the 46-kb TP53 deletion (below ASCAT's segment resolution, but independently 0 copies on per-exon depth) and the KDM3A amplicon amplitude (16 vs 18).",
  "Consequences: report purity as 0.57–0.67 and ploidy as 2.8–3.4; WGD, LOH landscape and gene-level states are [solid]; multiplicity-based claims that depend on the exact purity (e.g. RB1 p.Thr502Ile at 1.5–1.8 copies) are 'consistent with' rather than definitive. A third caller (Sequenza/PureCN, EC2) would break the tie on ploidy.")
writeLines(rep, file.path(cn, "ascat_report.md")); cat(rep, sep = "\n")
