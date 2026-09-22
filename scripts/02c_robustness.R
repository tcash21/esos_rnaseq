#!/usr/bin/env Rscript
# =============================================================================
# 02c_robustness.R — sensitivity analyses for the transcriptome-identity result (audit v0.1, items 14–16).
#   PROJ=$PWD Rscript scripts/02c_robustness.R          (uses the GDC-recipe quantification, 02b)
#
# 1. Similarity WITHOUT quantile mapping (raw log2 TPM) vs WITH — Spearman is rank-based so should be identical;
#    PCA projection is shown both ways.
# 2. Reference-only classifiers (nearest centroid, kNN) trained on cohort samples only, leave-one-out CV accuracy
#    for OS vs STS and per subtype, THEN the ESOS sample classified.
# 3. Bootstrap of gene selection: 500 random draws of 1,000 genes from the 5,000 most variable; fraction of draws
#    in which TARGET-OS is the top group by median rho / centroid rho.
# 4. Anatomical-site and age controls: cohort site/age distributions; posterior HOX genes by site
#    (TARGET-OS lower-limb vs pelvic/axial; TCGA-SARC retroperitoneum vs lower limb) vs the tumor.
# Output: results/02_expression/gdc_requant/robustness_report.md, pca_no_mapping.png, hox_by_site.png
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
proj <- Sys.getenv("PROJ", getwd()); xena <- file.path(proj, "refs/xena"); out <- file.path(proj, "results/02_expression/gdc_requant")
set.seed(42)

tum <- fread(file.path(proj, "results/02_expression/requant/ISM563041-2.gdc.genes.tsv"))[, .(gene = gene_name, tpm = tpm_unstranded)][, .(tpm = sum(tpm)), by = gene]
pm <- fread(file.path(xena, "gencode.v36.annotation.gtf.gene.probemap"))[, .(id, gene)]
load_cohort <- function(f) { m <- fread(cmd = paste("gzcat", shQuote(file.path(xena, f)))); setnames(m, 1, "id")
  m <- merge(pm, m, by = "id")[, id := NULL][, lapply(.SD, max), by = gene]
  keep <- grep("-01[AB]$", names(m), value = TRUE); mat <- as.matrix(m[, ..keep]); rownames(mat) <- m$gene; mat }
os <- load_cohort("TARGET-OS.star_tpm.tsv.gz"); sa <- load_cohort("TCGA-SARC.star_tpm.tsv.gz")
rd_clin <- function(f) fread(cmd = paste("gzcat", shQuote(file.path(xena, f))))
cs <- rd_clin("TCGA-SARC.clinical.tsv.gz"); co <- rd_clin("TARGET-OS.clinical.tsv.gz")
site_sa <- cs[match(colnames(sa), sample), tissue_or_organ_of_origin.diagnoses]
site_os <- co[match(colnames(os), sample), tissue_or_organ_of_origin.diagnoses]
age_sa <- cs[match(colnames(sa), sample), age_at_index.demographic]; age_os <- co[match(colnames(os), sample), age_at_index.demographic]
dx <- cs[match(colnames(sa), sample), primary_diagnosis.diagnoses]
sub_sa <- fcase(grepl("Leiomyosarcoma", dx), "LMS", grepl("Dedifferentiated liposarcoma", dx), "DDLPS",
                grepl("Undifferentiated|Malignant fibrous|Giant cell", dx), "UPS", grepl("Fibromyxosarcoma", dx), "MFS",
                grepl("nerve sheath", dx), "MPNST", grepl("Synovial", dx), "SS", default = "STS_other")

genes <- Reduce(intersect, list(tum$gene, rownames(os), rownames(sa)))
X <- cbind(os[genes, ], sa[genes, ]); X <- X[rowMeans(X) >= 1, ]; genes <- rownames(X)
grp <- c(rep("TARGET-OS", ncol(os)), sub_sa); coarse <- ifelse(grp == "TARGET-OS", "OS", "STS")
t_raw <- log2(tum[match(genes, gene), tpm] + 1); names(t_raw) <- genes
ref_q <- sort(rowMeans(apply(X, 2, sort))); t_q <- ref_q[rank(t_raw, ties.method = "average")]; names(t_q) <- genes
v <- apply(X, 1, var); top <- names(sort(v, decreasing = TRUE))[1:2000]

# ---------- 1. with vs without quantile mapping ----------
rho_of <- function(tv, g = top) apply(X[g, ], 2, function(s) cor(s, tv[g], method = "spearman"))
med_by <- function(r) sort(tapply(r, grp, median), decreasing = TRUE)
m_raw <- med_by(rho_of(t_raw)); m_q <- med_by(rho_of(t_q))
Xc <- t(X[top, ]); mu <- colMeans(Xc); pca <- prcomp(scale(Xc, center = mu, scale = FALSE), rank. = 5)
proj_pt <- function(tv) as.numeric((tv[top] - mu) %*% pca$rotation)
p_raw <- proj_pt(t_raw); p_q <- proj_pt(t_q)
pcs <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = grp)
pl <- ggplot(pcs, aes(PC1, PC2, colour = group)) + geom_point(alpha = .5, size = 1.6) +
  annotate("point", x = p_raw[1], y = p_raw[2], shape = 8, size = 6, stroke = 1.4) + annotate("label", x = p_raw[1], y = p_raw[2], label = "raw log2 TPM", vjust = -0.7, size = 3) +
  annotate("point", x = p_q[1], y = p_q[2], shape = 4, size = 6, stroke = 1.4) + annotate("label", x = p_q[1], y = p_q[2], label = "quantile-mapped", vjust = 1.7, size = 3) +
  theme_bw() + labs(title = "PCA projection of the ESOS tumor with and without quantile mapping")
ggsave(file.path(out, "pca_no_mapping.png"), pl, width = 9, height = 6, dpi = 130)
dist_to <- function(p) sort(sapply(split(seq_along(grp), grp), function(i) sqrt(sum((colMeans(pca$x[i, 1:5, drop = FALSE]) - p[1:5])^2))))

# ---------- 2. reference-only classifiers with LOO-CV ----------
Z <- X[top, ]
loo_centroid <- function(labels) { pred <- character(ncol(Z))
  for (j in seq_len(ncol(Z))) { cen <- sapply(split(setdiff(seq_len(ncol(Z)), j), labels[-j]), function(i) rowMeans(Z[, i, drop = FALSE]))
    pred[j] <- names(which.max(apply(cen, 2, function(c) cor(c, Z[, j], method = "spearman")))) }
  pred }
knn_pred <- function(labels, k = 15) { R <- cor(Z, method = "spearman"); diag(R) <- NA
  sapply(seq_len(ncol(Z)), function(j) { nn <- order(R[, j], decreasing = TRUE, na.last = NA)[1:k]; names(which.max(table(labels[nn]))) }) }
acc <- function(pred, labels) mean(pred == labels)
cv <- list(centroid_coarse = acc(loo_centroid(coarse), coarse), knn_coarse = acc(knn_pred(coarse), coarse),
           centroid_subtype = acc(loo_centroid(grp), grp), knn_subtype = acc(knn_pred(grp), grp))
# classify the ESOS sample (trained on all reference samples)
cen_all <- sapply(split(seq_len(ncol(Z)), grp), function(i) rowMeans(Z[, i, drop = FALSE]))
esos_centroid <- names(which.max(apply(cen_all, 2, function(c) cor(c, t_raw[top], method = "spearman"))))
rr <- rho_of(t_raw); esos_knn <- names(which.max(table(grp[order(rr, decreasing = TRUE)[1:15]])))
knn_votes <- table(grp[order(rr, decreasing = TRUE)[1:15]])

# ---------- 3. bootstrap gene selection ----------
top5k <- names(sort(v, decreasing = TRUE))[1:5000]
boot <- replicate(500, { g <- sample(top5k, 1000); r <- rho_of(t_raw, g)
  c(median_top = names(which.max(tapply(r, grp, median))),
    centroid_top = names(which.max(apply(sapply(split(seq_along(grp), grp), function(i) rowMeans(X[g, i, drop = FALSE])), 2, function(c) cor(c, t_raw[g], method = "spearman")))),
    gap = unname(sort(tapply(r, grp, median), decreasing = TRUE)[1] - sort(tapply(r, grp, median), decreasing = TRUE)[2])) })
boot_med <- mean(boot["median_top", ] == "TARGET-OS"); boot_cen <- mean(boot["centroid_top", ] == "TARGET-OS")
gap_q <- quantile(as.numeric(boot["gap", ]), c(.025, .5, .975))

# ---------- 4. site / age controls and HOX by site ----------
site_grp <- c(fcase(grepl("lower limb", site_os), "OS lower limb", grepl("upper limb", site_os), "OS upper limb",
                    grepl("Pelvic|skull|Bone, NOS|Short bones", site_os), "OS pelvic/axial/other", default = "OS site n/r"),
              fcase(grepl("Retroperitoneum", site_sa), "STS retroperitoneum", grepl("lower limb", site_sa), "STS lower limb",
                    grepl("Uterus|Myometrium", site_sa), "STS uterus", default = "STS other/trunk"))
hox <- c("HOXA9", "HOXA10", "HOXA11", "HOXA13", "HOXC10", "HOXC11", "HOXD10", "HOXD11", "PITX1", "TBX5", "MSX1", "SATB2", "SP7", "RUNX2")
hox <- hox[hox %in% genes]
hl <- rbindlist(lapply(hox, function(g) data.table(gene = g, site = site_grp, expr = X[g, ])))
hl_t <- data.table(gene = hox, expr = t_raw[hox])
hs <- hl[, .(n = .N, median = round(median(expr), 2), q25 = round(quantile(expr, .25), 2), q75 = round(quantile(expr, .75), 2)), by = .(gene, site)]
hs <- merge(hs, hl_t[, .(gene, tumor = round(expr, 2))], by = "gene")
hs[, pct := { g0 <- gene; s0 <- site; t0 <- tumor
  round(100 * sapply(seq_len(.N), function(i) mean(hl[gene == g0[i] & site == s0[i], expr] <= t0[i]))) }]
# full HOX clusters vs TARGET-OS (is the reduction cluster-specific?)
hoxall <- grep("^HOX[ABCD][0-9]+$", genes, value = TRUE)
i_os <- grp == "TARGET-OS"; i_rp <- site_grp == "STS retroperitoneum"
hz <- data.table(gene = hoxall, cluster = substr(hoxall, 4, 4), paralog = as.integer(sub("HOX[ABCD]", "", hoxall)),
                 tumor = round(t_raw[hoxall], 2), OS_median = round(apply(X[hoxall, i_os], 1, median), 2),
                 z_OS = round((t_raw[hoxall] - rowMeans(X[hoxall, i_os])) / apply(X[hoxall, i_os], 1, sd), 2),
                 RP_STS_median = round(apply(X[hoxall, i_rp], 1, median), 2),
                 z_RP_STS = round((t_raw[hoxall] - rowMeans(X[hoxall, i_rp])) / apply(X[hoxall, i_rp], 1, sd), 2))[order(cluster, paralog)]
ph <- ggplot(hl[gene %in% c("HOXA9", "HOXA10", "HOXA11", "HOXC10", "PITX1", "SATB2", "SP7")], aes(site, expr, fill = grepl("^OS", site))) +
  geom_boxplot(outlier.size = .4, show.legend = FALSE) + geom_hline(data = hl_t[gene %in% c("HOXA9", "HOXA10", "HOXA11", "HOXC10", "PITX1", "SATB2", "SP7")], aes(yintercept = expr), colour = "red") +
  facet_wrap(~gene, scales = "free_y", ncol = 4) + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "log2(TPM+1)", title = "Posterior HOX / osteoblastic genes by cohort and anatomical site; red = ESOS tumor")
ggsave(file.path(out, "hox_by_site.png"), ph, width = 12, height = 6, dpi = 130)

# ---------- report ----------
f2 <- function(x) sprintf("%.2f", x)
rep <- c("# Robustness of the transcriptome-identity result", paste0("_generated ", Sys.Date(), "; GDC-recipe quantification; audit items 14–16_"), "",
  "## 1. Quantile mapping does not drive the result",
  sprintf("- Median Spearman rho by group, raw log2 TPM: %s", paste(names(m_raw), f2(m_raw), collapse = ", ")),
  sprintf("- Median Spearman rho by group, quantile-mapped: %s", paste(names(m_q), f2(m_q), collapse = ", ")),
  "  (identical, as expected: Spearman is invariant to the monotone mapping)",
  sprintf("- PCA distance to group centroid (PC1–5), raw: %s", paste(names(dist_to(p_raw)), f2(dist_to(p_raw)), collapse = ", ")),
  sprintf("- PCA distance to group centroid (PC1–5), mapped: %s", paste(names(dist_to(p_q)), f2(dist_to(p_q)), collapse = ", ")),
  "  ![PCA](pca_no_mapping.png)", "",
  "## 2. Reference-only classifiers (trained on cohort samples only; leave-one-out CV)",
  sprintf("- Nearest-centroid, OS vs STS: LOO accuracy %.1f%%; by subtype: %.1f%%", 100 * cv$centroid_coarse, 100 * cv$centroid_subtype),
  sprintf("- kNN (k=15), OS vs STS: LOO accuracy %.1f%%; by subtype: %.1f%%", 100 * cv$knn_coarse, 100 * cv$knn_subtype),
  sprintf("- ESOS sample, classified after training: nearest-centroid → **%s**; kNN → **%s** (votes: %s)", esos_centroid, esos_knn, paste(names(knn_votes), knn_votes, collapse = ", ")), "",
  "## 3. Bootstrap of gene selection (500 draws of 1,000 genes from the 5,000 most variable)",
  sprintf("- TARGET-OS is the top group by median rho in **%.1f%%** of draws and by centroid rho in **%.1f%%**", 100 * boot_med, 100 * boot_cen),
  sprintf("- Margin between the top group and the runner-up (median rho): median %s, 95%% interval %s–%s", f2(gap_q[2]), f2(gap_q[1]), f2(gap_q[3])), "",
  "## 4. Cohort composition (confounders the audit raised)",
  sprintf("- TARGET-OS: n=%d, age median %s y (IQR %s–%s); sites: %s", ncol(os), median(age_os, na.rm = TRUE), quantile(age_os, .25, na.rm = TRUE), quantile(age_os, .75, na.rm = TRUE),
          paste(names(table(site_grp[grp == "TARGET-OS"])), table(site_grp[grp == "TARGET-OS"]), collapse = "; ")),
  sprintf("- TCGA-SARC: n=%d, age median %s y (IQR %s–%s); sites: %s", ncol(sa), median(age_sa, na.rm = TRUE), quantile(age_sa, .25, na.rm = TRUE), quantile(age_sa, .75, na.rm = TRUE),
          paste(names(table(site_grp[grp != "TARGET-OS"])), table(site_grp[grp != "TARGET-OS"]), collapse = "; ")),
  "- The ESOS patient (mid-60s, retroperitoneum) is closer in age and site to TCGA-SARC than to TARGET-OS; the age/site confound therefore works *against* an osteosarcoma call, not for it.", "",
  "## 5. Posterior HOX and osteoblastic genes by anatomical site (log2 TPM+1; %ile = fraction of that site group at or below the tumor)",
  "| gene | site group | n | median | IQR | tumor | %ile |", "|---|---|---|---|---|---|---|",
  hs[order(gene, site), sprintf("| %s | %s | %d | %s | %s–%s | %s | %d |", gene, site, n, f2(median), f2(q25), f2(q75), f2(tumor), pct)], "",
  "  ![HOX](hox_by_site.png)", "",
  "## 6. All four HOX clusters: tumor vs TARGET-OS and vs retroperitoneal TCGA-SARC (z-scores)",
  "| gene | tumor | OS median | z vs OS | RP-STS median | z vs RP-STS |", "|---|---|---|---|---|---|",
  hz[, sprintf("| %s | %s | %s | %s | %s | %s |", gene, f2(tumor), f2(OS_median), f2(z_OS), f2(RP_STS_median), f2(z_RP_STS))], "",
  "## Reading guide",
  "- If HOXA10/HOXC10/PITX1 are low in the tumor relative to *both* lower-limb OS and retroperitoneal STS, the reduction is not explained by anatomical site; if retroperitoneal STS are also low, it is a site effect.",
  "- Pelvic/axial OS in TARGET is a small group (n≈6); treat that comparison as indicative only.")
writeLines(rep, file.path(out, "robustness_report.md")); cat(rep, sep = "\n")
