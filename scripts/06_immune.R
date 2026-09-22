#!/usr/bin/env Rscript
# =============================================================================
# 06_immune.R — immune / stromal composition of the tumor vs TARGET-OS and TCGA-SARC (Mac).
#   Rscript scripts/06_immune.R
#
# Two methods on the SAME quantification (GDC recipe, 02b) as the cohorts:
#   MCP-counter  rank-based marker scores; only comparable across samples -> run on tumor + cohorts together
#   quanTIseq    constrained regression on TPM (TIL10 signature) -> absolute cell fractions per sample
# Inputs : results/02_expression/requant/ISM563041-2.gdc.genes.tsv, refs/xena/{TARGET-OS,TCGA-SARC}.star_tpm.tsv.gz
# Outputs: results/06_hla/immune_report.md, immune_scores.tsv, immune_mcp.png
# Caveats: purity 0.57 (FACETS) means ~43% of RNA is non-tumor by construction; cohort purities vary too.
#          Tumor is total RNA, cohorts poly-A (library prep); quanTIseq tolerates this less well than MCP-counter.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(MCPcounter); library(quantiseqr) })
proj <- Sys.getenv("PROJ", getwd()); xena <- file.path(proj, "refs/xena"); out <- file.path(proj, "results/06_hla")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---------- expression: tumor (GDC TPM) + cohorts (log2(TPM+1) -> TPM) on gene symbols ----------
tum <- fread(file.path(proj, "results/02_expression/requant/ISM563041-2.gdc.genes.tsv"))[gene_type == "protein_coding", .(tpm = sum(tpm_unstranded)), by = .(gene = gene_name)]
pm <- fread(file.path(xena, "gencode.v36.annotation.gtf.gene.probemap"))[, .(id, gene)]
load_cohort <- function(f) {
  m <- fread(cmd = paste("gzcat", shQuote(file.path(xena, f)))); setnames(m, 1, "id")
  m <- merge(pm, m, by = "id")[, id := NULL][, lapply(.SD, max), by = gene]
  keep <- grep("-01[AB]$", names(m), value = TRUE); mat <- 2^as.matrix(m[, ..keep]) - 1; rownames(mat) <- m$gene; mat
}
os <- load_cohort("TARGET-OS.star_tpm.tsv.gz"); sa <- load_cohort("TCGA-SARC.star_tpm.tsv.gz")
genes <- Reduce(intersect, list(tum$gene, rownames(os), rownames(sa)))
X <- cbind(ESOS = tum[match(genes, gene), tpm], os[genes, ], sa[genes, ])
X <- sweep(X, 2, colSums(X), "/") * 1e6                      # re-normalise to TPM over the shared gene space
grp <- c("ESOS", rep("TARGET-OS", ncol(os)), rep("TCGA-SARC", ncol(sa)))
cat(sprintf("shared protein-coding genes: %d; samples: %d\n", length(genes), ncol(X)))

# ---------- MCP-counter (log2 TPM, all samples together) ----------
mcp <- MCPcounter.estimate(log2(X + 1), featuresType = "HUGO_symbols")
mcp_dt <- as.data.table(t(mcp), keep.rownames = "sample")[, group := grp]
mcp_long <- melt(mcp_dt, id.vars = c("sample", "group"), variable.name = "population", value.name = "score")
mcp_sum <- mcp_long[group != "ESOS", .(os_pct = if (group[1] == "TARGET-OS") NA else NA), by = population]   # placeholder, replaced below
pct_in <- function(pop, g) { v <- mcp_long[population == pop & group == g, score]; t <- mcp_long[population == pop & group == "ESOS", score]; round(100 * mean(v <= t)) }
mcp_tab <- mcp_long[group == "ESOS", .(population, tumor_score = round(score, 2))][
  , `:=`(pct_TARGET_OS = sapply(population, pct_in, g = "TARGET-OS"), pct_TCGA_SARC = sapply(population, pct_in, g = "TCGA-SARC"))]

p <- ggplot(mcp_long[group != "ESOS"], aes(group, score, fill = group)) + geom_boxplot(outlier.size = .4, show.legend = FALSE) +
  geom_hline(data = mcp_long[group == "ESOS"], aes(yintercept = score), colour = "red", linewidth = .9) +
  facet_wrap(~population, scales = "free_y", ncol = 5) + theme_bw(base_size = 9) +
  labs(x = NULL, y = "MCP-counter score (log2)", title = "Immune / stromal populations: cohorts (boxes) vs ESOS tumor (red line)")
ggsave(file.path(out, "immune_mcp.png"), p, width = 11, height = 5.5, dpi = 130)

# ---------- quanTIseq (TPM, per-sample absolute fractions) ----------
qt <- run_quantiseq(X, signature_matrix = "TIL10", is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE)
qt_dt <- as.data.table(qt)[, group := grp]
qt_long <- melt(qt_dt, id.vars = c("Sample", "group"), variable.name = "cell", value.name = "fraction")
qt_tab <- qt_long[group == "ESOS", .(cell, tumor_fraction = round(100 * fraction, 1))][
  , `:=`(TARGET_OS_median = sapply(cell, function(c) round(100 * median(qt_long[cell == c & group == "TARGET-OS", fraction]), 1)),
         pct_TARGET_OS = sapply(cell, function(c) round(100 * mean(qt_long[cell == c & group == "TARGET-OS", fraction] <= qt_long[cell == c & group == "ESOS", fraction]))),
         TCGA_SARC_median = sapply(cell, function(c) round(100 * median(qt_long[cell == c & group == "TCGA-SARC", fraction]), 1)))]

fwrite(rbind(mcp_tab[, .(method = "MCP-counter", feature = population, tumor = tumor_score, pct_TARGET_OS, pct_TCGA_SARC)],
             qt_tab[, .(method = "quanTIseq", feature = cell, tumor = tumor_fraction, pct_TARGET_OS, pct_TCGA_SARC = NA)]),
       file.path(out, "immune_scores.tsv"), sep = "\t")

# ---------- report ----------
rep <- c("# Immune / stromal composition — ESOS tumor vs TARGET-OS and TCGA-SARC",
  paste0("_generated ", Sys.Date(), "; GDC-recipe quantification for all samples; research-grade_"), "",
  "## MCP-counter (rank scores; percentile = fraction of cohort samples at or below the tumor)",
  "| population | tumor score | %ile in TARGET-OS | %ile in TCGA-SARC |", "|---|---|---|---|",
  mcp_tab[, sprintf("| %s | %.2f | %d | %d |", population, tumor_score, pct_TARGET_OS, pct_TCGA_SARC)], "",
  "![MCP](immune_mcp.png)", "",
  "## quanTIseq (TIL10; estimated % of cells)",
  "| cell type | tumor % | TARGET-OS median % | %ile in TARGET-OS | TCGA-SARC median % |", "|---|---|---|---|---|",
  qt_tab[, sprintf("| %s | %.1f | %.1f | %d | %.1f |", cell, tumor_fraction, TARGET_OS_median, pct_TARGET_OS, TCGA_SARC_median)], "",
  "## Reading guide",
  "- FACETS purity 0.57 -> ~43% non-tumor cells; quanTIseq's 'Other' should be read against that (tumor + stroma + unassigned).",
  "- Consistent with the marker-level view (step 02): macrophage-rich, T-cell-poor. Whether the fibroblast/endothelial signal is abdominal stroma or tumor is the open ESOS question.",
  "- Library prep differs (total RNA vs poly-A); MCP-counter percentiles are more robust to this than quanTIseq fractions.")
writeLines(rep, file.path(out, "immune_report.md")); cat(rep, sep = "\n")
