#!/usr/bin/env Rscript
# facets_fit.R <pileup.csv.gz> <outdir> <cval> [<cval> ...]
# Runs FACETS (Shen & Seshan 2016) once per cval and writes, for each:
#   facets_cval<C>_segments.tsv  chrom, start, end, num.mark, cnlr.median, mafR, tcn.em, lcn.em, cf.em
#   facets_cval<C>_fit.tsv       purity, ploidy, dipLogR, loglik, flags
#   facets_cval<C>.png           standard FACETS diagnostic plot
# The Mac-side 04_copynumber.R consumes these (WGD call, LOH fraction, HRD scores, gene-level table).
suppressPackageStartupMessages(library(facets))
args <- commandArgs(trailingOnly = TRUE)
pileup <- args[1]; out <- args[2]; cvals <- as.integer(args[-(1:2)])
set.seed(1234)

rc <- readSnpMatrix(pileup)
cat(sprintf("pileup: %d SNP positions read\n", nrow(rc)))

for (cv in cvals) {
  # ndepth 25: exome depth is high, so demand real coverage; snp.nbhd 250 to thin dense exon SNPs
  pre <- preProcSample(rc, gbuild = "hg19", ndepth = 25, snp.nbhd = 250, cval = 25)
  oo  <- procSample(pre, cval = cv)
  fit <- emcncf(oo)
  seg <- fit$cncf
  seg$chrom <- ifelse(seg$chrom == 23, "X", as.character(seg$chrom))
  keep <- c("chrom", "start", "end", "num.mark", "nhet", "cnlr.median", "mafR", "tcn.em", "lcn.em", "cf.em")
  write.table(seg[, intersect(keep, names(seg))], file.path(out, sprintf("facets_cval%d_segments.tsv", cv)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(cval = cv, purity = fit$purity, ploidy = fit$ploidy, dipLogR = fit$dipLogR,
                         loglik = fit$loglik, flags = paste(oo$flags, collapse = ";"), n_segments = nrow(seg)),
              file.path(out, sprintf("facets_cval%d_fit.tsv", cv)), sep = "\t", quote = FALSE, row.names = FALSE)
  png(file.path(out, sprintf("facets_cval%d.png", cv)), width = 1600, height = 1000, res = 120)
  plotSample(x = oo, emfit = fit, sname = sprintf("ESOS tumor  cval=%d  purity=%.2f  ploidy=%.2f", cv, fit$purity, fit$ploidy))
  dev.off()
  cat(sprintf("cval %d: purity %.3f  ploidy %.2f  dipLogR %.3f  segments %d  flags: %s\n",
              cv, fit$purity, fit$ploidy, fit$dipLogR, nrow(seg), paste(oo$flags, collapse = ";")))
}
