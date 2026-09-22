#!/usr/bin/env Rscript
# =============================================================================
# 02_expression.R — does the tumor's transcriptome look like bone osteosarcoma or like a
# soft-tissue sarcoma?  Runs on the Mac.      Rscript scripts/02_expression.R
#
# Inputs : data/sema4/**/ISM563041-2.genes.results     (RSEM, RefSeq, TPM)  — the tumor
#          refs/xena/TARGET-OS.star_tpm.tsv.gz         (88 pediatric OS, log2(TPM+1), GDC STAR)
#          refs/xena/TCGA-SARC.star_tpm.tsv.gz         (258 soft-tissue sarcomas, same pipeline)
#          refs/xena/TCGA-SARC.clinical.tsv.gz         (histology labels)
#          refs/xena/gencode.v36.annotation.gtf.gene.probemap (ENSG -> symbol)
# Pipeline caveat: the tumor is RSEM/RefSeq, the cohorts are STAR/GENCODE. Everything below is
# rank-based: the tumor is quantile-mapped onto the pooled cohort distribution, and
# comparisons use Spearman correlation, centroid ranks and per-gene z-scores of the mapped values.
# A residual pipeline effect cannot be excluded; it would push the tumor away from BOTH cohorts
# equally, not toward one of them, so the OS-vs-STS *direction* is robust to it.
#
# Outputs: results/02_expression/
#   expression_report.md, pca_cohorts.png, correlation_by_group.png,
#   gene_zscores.tsv (key genes vs TARGET-OS and vs TCGA-SARC), tumor_extremes_vs_TARGETOS.tsv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
proj <- Sys.getenv("PROJ", getwd())
xena <- file.path(proj, "refs/xena")
# QUANT=sema4 (default): Sema4's RSEM/RefSeq TPM.  QUANT=gdc: our GDC-recipe re-quantification (02b), same
# STAR/GENCODE v36 pipeline as the cohorts -> like-for-like (library prep still differs: total RNA vs poly-A).
quant <- Sys.getenv("QUANT", "sema4")
out <- file.path(proj, "results/02_expression", if (quant == "gdc") "gdc_requant" else "."); dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# ---------- tumor ----------
if (quant == "gdc") {
  gf <- file.path(proj, "results/02_expression/requant/ISM563041-2.gdc.genes.tsv")
  tum <- fread(gf)[, .(gene = gene_name, tpm = tpm_unstranded)][, .(tpm = sum(tpm)), by = gene]   # unstranded = GDC's cohort convention
} else {
  gf <- list.files(file.path(proj, "data/sema4"), "ISM563041-2\\.genes\\.results$", recursive = TRUE, full.names = TRUE)[1]
  tum <- fread(gf)[, .(gene = sub("^HGNC:[0-9]+_", "", gene_id), tpm = TPM)][, .(tpm = sum(tpm)), by = gene]
}
cat(sprintf("tumor quantification: %s (%s)\n", quant, gf))

# ---------- cohorts ----------
pm <- fread(file.path(xena, "gencode.v36.annotation.gtf.gene.probemap"))[, .(id, gene)]
load_cohort <- function(f) {
  m <- fread(cmd = paste("gzcat", shQuote(file.path(xena, f))))
  setnames(m, 1, "id"); m <- merge(pm, m, by = "id")[, id := NULL]
  m <- m[, lapply(.SD, max), by = gene]                      # collapse duplicate symbols
  keep <- grep("-01[AB]$", names(m), value = TRUE)          # primary tumours only
  list(genes = m$gene, mat = as.matrix(m[, ..keep]))
}
os <- load_cohort("TARGET-OS.star_tpm.tsv.gz"); sa <- load_cohort("TCGA-SARC.star_tpm.tsv.gz")
rownames(os$mat) <- os$genes; rownames(sa$mat) <- sa$genes

# TCGA-SARC histology
cl <- fread(cmd = paste("gzcat", shQuote(file.path(xena, "TCGA-SARC.clinical.tsv.gz"))), select = c("sample", "primary_diagnosis.diagnoses"))
setnames(cl, c("sample", "dx"))
cl[, group := fcase(grepl("Leiomyosarcoma", dx), "LMS", grepl("Dedifferentiated liposarcoma", dx), "DDLPS",
                    grepl("Undifferentiated|Malignant fibrous histiocytoma|Giant cell", dx), "UPS",
                    grepl("Fibromyxosarcoma", dx), "MFS", grepl("nerve sheath", dx), "MPNST",
                    grepl("Synovial", dx), "SS", default = "STS_other")]
sa_group <- cl[match(colnames(sa$mat), sample), group]; sa_group[is.na(sa_group)] <- "STS_other"

# ---------- common gene space ----------
genes <- Reduce(intersect, list(tum$gene, rownames(os$mat), rownames(sa$mat)))
X <- cbind(os$mat[genes, ], sa$mat[genes, ])                # cohorts, log2(TPM+1)
grp <- c(rep("TARGET-OS", ncol(os$mat)), sa_group)
expressed <- rowMeans(X) >= 1                               # drop genes silent in the cohorts
X <- X[expressed, ]; genes <- rownames(X)
t_log <- log2(tum[match(genes, gene), tpm] + 1); names(t_log) <- genes

# quantile-map the tumor onto the pooled cohort distribution (rank -> reference quantile)
ref_q <- sort(rowMeans(apply(X, 2, sort)))                  # mean quantile profile of the cohorts
t_q <- ref_q[rank(t_log, ties.method = "average")]; names(t_q) <- genes
cat(sprintf("genes in common & expressed: %d; TARGET-OS n=%d; TCGA-SARC n=%d (%s)\n", length(genes), ncol(os$mat), ncol(sa$mat),
            paste(names(table(sa_group)), table(sa_group), sep = "=", collapse = ", ")))

# ---------- 1. correlation with every cohort sample, top-variable genes ----------
v <- apply(X, 1, var); top <- names(sort(v, decreasing = TRUE))[1:2000]
rho <- apply(X[top, ], 2, function(s) cor(s, t_q[top], method = "spearman"))
cor_dt <- data.table(sample = colnames(X), group = grp, rho = rho)
cor_sum <- cor_dt[, .(n = .N, median_rho = median(rho), q25 = quantile(rho, .25), q75 = quantile(rho, .75)), by = group][order(-median_rho)]
nn <- cor_dt[order(-rho)][1:25, .N, by = group][order(-N)]

# ---------- 2. centroid classification ----------
cent <- sapply(split(seq_along(grp), grp), function(i) rowMeans(X[top, i, drop = FALSE]))
cent_rho <- sort(apply(cent, 2, function(c) cor(c, t_q[top], method = "spearman")), decreasing = TRUE)

# ---------- 3. PCA of cohorts, tumor projected ----------
Xc <- t(X[top, ]); mu <- colMeans(Xc); pca <- prcomp(scale(Xc, center = mu, scale = FALSE), rank. = 5)
t_pc <- (t_q[top] - mu) %*% pca$rotation
pcs <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = grp)
ve <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2))
p <- ggplot(pcs, aes(PC1, PC2, colour = group)) + geom_point(alpha = .6, size = 1.8) +
  annotate("point", x = t_pc[1], y = t_pc[2], shape = 8, size = 6, stroke = 1.5, colour = "black") +
  annotate("label", x = t_pc[1], y = t_pc[2], label = "ESOS tumor", vjust = -0.8, size = 3.5) +
  labs(x = sprintf("PC1 (%d%%)", ve[1]), y = sprintf("PC2 (%d%%)", ve[2]),
       title = "TARGET-OS + TCGA-SARC, top 2000 variable genes; tumor quantile-mapped and projected") + theme_bw()
ggsave(file.path(out, "pca_cohorts.png"), p, width = 9, height = 6, dpi = 130)
# distance in PC space to each group's centroid (Mahalanobis-lite: euclidean on PC1-5)
pc5 <- pca$x[, 1:5]; t5 <- as.numeric(t_pc[1:5])
pc_dist <- sort(sapply(split(seq_along(grp), grp), function(i) sqrt(sum((colMeans(pc5[i, , drop = FALSE]) - t5)^2))))

p2 <- ggplot(cor_dt, aes(reorder(group, rho, median), rho, fill = group)) + geom_boxplot(outlier.size = .6, show.legend = FALSE) +
  coord_flip() + labs(x = NULL, y = "Spearman rho with the ESOS tumor (top 2000 variable genes)") + theme_bw()
ggsave(file.path(out, "correlation_by_group.png"), p2, width = 7, height = 4, dpi = 130)

# ---------- 4. gene-level z-scores ----------
key <- list(
  osteoblastic = c("RUNX2", "SP7", "SPP1", "IBSP", "ALPL", "COL1A1", "BGLAP", "SATB2", "DLX5", "MEPE", "DMP1"),
  cell_cycle   = c("CDK4", "CDK6", "CCND1", "CCND3", "CCNE1", "CDKN2A", "CDKN1A", "E2F1", "E2F3", "MYC", "MKI67", "TOP2A", "AURKB", "COPS3"),
  amplified    = c("KDM3A", "MAP4K4"),
  tp53_rb      = c("TP53", "RB1", "MDM2", "CDKN1C", "PTEN"),
  telomere_alt = c("TERT", "ATRX", "DAXX"),
  immune       = c("PTPRC", "CD3E", "CD8A", "CD4", "FOXP3", "NKG7", "CD274", "PDCD1", "CTLA4", "CD68", "CD163", "CD276", "B2M", "HLA-A", "HLA-B", "TAP1"),
  sts_markers  = c("DES", "ACTA2", "MYOCD", "CNN1", "S100B", "SOX10", "TLE1", "MDM2", "HMGA2", "PLAG1"),
  rtk_targets  = c("IGF1R", "PDGFRA", "PDGFRB", "KIT", "MET", "VEGFA", "FGFR1", "ERBB2", "GD2_B4GALNT1" = "B4GALNT1"))
kg <- unique(unlist(key)); kg <- kg[kg %in% genes]
z <- function(g, idx) (t_q[g] - rowMeans(X[g, idx, drop = FALSE])) / apply(X[g, idx, drop = FALSE], 1, sd)
pct <- function(g, idx) sapply(g, function(x) mean(X[x, idx] <= t_q[x]))
i_os <- grp == "TARGET-OS"; i_sa <- grp != "TARGET-OS"
gz <- data.table(gene = kg, set = sapply(kg, function(g) names(key)[sapply(key, function(s) g %in% s)][1]),
                 tumor_tpm = tum[match(kg, gene), tpm], mapped_log2 = round(t_q[kg], 2),
                 z_vs_OS = round(z(kg, i_os), 2), pct_vs_OS = round(100 * pct(kg, i_os)),
                 z_vs_STS = round(z(kg, i_sa), 2), pct_vs_STS = round(100 * pct(kg, i_sa)),
                 OS_median_log2 = round(apply(X[kg, i_os], 1, median), 2), STS_median_log2 = round(apply(X[kg, i_sa], 1, median), 2))
fwrite(gz, file.path(out, "gene_zscores.tsv"), sep = "\t")

# genome-wide: where is the tumor extreme relative to TARGET-OS?
# Exclude quantification-sensitive classes (ribosomal proteins, miRNA/snoRNA hosts, antisense, pseudogenes) and genes with
# near-zero cohort variance, and require a >= 4-fold shift: otherwise pipeline differences dominate the list.
sd_os <- apply(X[, i_os], 1, sd)
zall <- (t_q - rowMeans(X[, i_os])) / sd_os
ext <- data.table(gene = genes, z_vs_OS = round(zall, 2), tumor_mapped = round(t_q, 2), OS_mean = round(rowMeans(X[, i_os]), 2), OS_sd = round(sd_os, 2),
                  z_vs_STS = round((t_q - rowMeans(X[, i_sa])) / apply(X[, i_sa], 1, sd), 2))
# Library-prep classes (tumor = total RNA, cohorts = poly-A): histone mRNAs, sn/sno/scaRNA, 7SL, Y_RNA, MT- genes.
# PAR1 genes: present on X and Y in the GENCODE reference -> multimapped -> dropped by STAR GeneCounts in our run.
par1 <- c("PLCXD1","GTPBP6","PPP2R3B","SHOX","CRLF2","CSF2RA","IL3RA","SLC25A6","ASMTL","P2RY8","AKAP17A","ASMT","DHRSX","ZBED1","CD99","XG")
ext <- ext[!grepl("^(RPL|RPS|MRPL|MRPS|MIR|SNOR|SNHG|SCARNA|LINC|LOC|RNU|RN7S|Y_RNA|U[0-9]$|MT-|H1-|H[1-4]C[0-9]|H2[AB]C[0-9]|H2[AB]Z|A[CLP][0-9]{6}\\.|Z[0-9]{5}\\.)|-AS[0-9]*$|P[0-9]+$", gene) &
           !gene %in% par1 & OS_sd >= 0.5 & abs(tumor_mapped - OS_mean) >= 2][order(-abs(z_vs_OS))]
fwrite(ext, file.path(out, "tumor_extremes_vs_TARGETOS.tsv"), sep = "\t")

# ---------- report ----------
f2 <- function(x) sprintf("%.2f", x)
gz_tab <- function(s) c("| gene | tumor TPM | z vs OS | %ile OS | z vs STS | %ile STS |", "|---|---|---|---|---|---|",
  gz[set == s][order(-z_vs_OS), sprintf("| %s | %.1f | %s | %d | %s | %d |", gene, tumor_tpm, f2(z_vs_OS), pct_vs_OS, f2(z_vs_STS), pct_vs_STS)])
rep <- c("# Expression — ESOS tumor vs TARGET-OS (bone OS) and TCGA-SARC (soft-tissue sarcoma)",
  paste0("_generated ", Sys.Date(), "; tumor quantification: **", if (quant == "gdc") "GDC recipe (STAR 2-pass, GRCh38, GENCODE v36; step 02b) — same pipeline as the cohorts" else "Sema4 RSEM/RefSeq (cross-pipeline; see caveat)", "**; rank-based_"), "",
  sprintf("Genes compared: %d. Cohorts: TARGET-OS n=%d; TCGA-SARC n=%d (%s).", length(genes), sum(i_os), sum(i_sa),
          paste(names(table(sa_group)), table(sa_group), sep = "=", collapse = ", ")), "",
  "## 1. Which cohort does the tumor resemble? (Spearman, top 2000 variable genes)",
  "| group | n | median rho | IQR |", "|---|---|---|---|",
  cor_sum[, sprintf("| %s | %d | **%s** | %s–%s |", group, n, f2(median_rho), f2(q25), f2(q75))], "",
  sprintf("- 25 nearest cohort samples by correlation: %s", nn[, paste0(group, " ", N, collapse = ", ")]),
  sprintf("- Centroid correlation, best to worst: %s", paste(names(cent_rho), f2(cent_rho), collapse = ", ")),
  sprintf("- PCA (PC1–5) distance to group centroid, nearest first: %s", paste(names(pc_dist), f2(pc_dist), collapse = ", ")), "",
  "![PCA](pca_cohorts.png)", "", "![correlation](correlation_by_group.png)", "",
  "## 2. Key genes (z-score and percentile of the tumor within each cohort)",
  "### Osteoblastic program", gz_tab("osteoblastic"), "", "### Cell cycle / OS amplicons", gz_tab("cell_cycle"), "",
  "### Focal amplification / kataegis genes", gz_tab("amplified"), "", "### TP53–RB axis", gz_tab("tp53_rb"), "",
  "### Telomere maintenance", gz_tab("telomere_alt"), "", "### Immune", gz_tab("immune"), "",
  "### Soft-tissue-sarcoma lineage markers (LMS: DES/ACTA2/MYOCD/CNN1; MPNST: S100B/SOX10; SS: TLE1; DDLPS: MDM2/HMGA2)", gz_tab("sts_markers"), "",
  "### Receptor tyrosine kinases / drug targets", gz_tab("rtk_targets"), "",
  if (quant == "gdc") c("## 3. Genome-wide extremes vs TARGET-OS (same quantification pipeline; library-prep-sensitive classes and PAR1 genes excluded; cohort SD >= 0.5; >= 4-fold)",
    "_Remaining confounder: total-RNA (tumor) vs poly-A (cohorts) library prep. Treat single genes as leads, not results._") else
  c("## 3. Genome-wide extremes vs TARGET-OS — NOT RELIABLE for this cross-pipeline run",
    "_Dominated by housekeeping/paralog genes = RSEM/RefSeq vs STAR/GENCODE differences. See results/02_expression/gdc_requant/ for the like-for-like run (QUANT=gdc)._"),
  "Up:", ext[z_vs_OS > 0][1:25, paste0(gene, " (", f2(z_vs_OS), ")", collapse = ", ")], "",
  "Down:", ext[z_vs_OS < 0][1:25, paste0(gene, " (", f2(z_vs_OS), ")", collapse = ", ")], "",
  "## Reading guide",
  "- If median rho and centroid/PCA distance all favour TARGET-OS over every STS subtype, the transcriptome is bone-OS-like; if UPS/MFS win, it is STS-like with osteoid.",
  "- z vs OS near 0 = typical osteosarcoma; |z| > 2 = unusual even for OS. Percentiles are the fraction of cohort samples at or below the tumor.",
  "- The tumor is 57% pure (FACETS); immune/stromal genes are diluted by design and cohort samples vary in purity too.")
writeLines(rep, file.path(out, "expression_report.md")); cat(rep, sep = "\n")
