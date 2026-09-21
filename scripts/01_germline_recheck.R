#!/usr/bin/env Rscript
# =============================================================================
# 01_germline_recheck.R — triage the re-annotated germline exome.
# Run after 01_germline_recheck.sh:   Rscript 01_germline_recheck.R
#
# Produces results/01_germline/
#   tier1_pathogenic_exomewide.tsv   ClinVar P/LP anywhere in the exome (>=1 star)
#   tier2_panel_rare.tsv             rare variants in cancer-predisposition genes
#                                    (protein-altering if a consequence annotation exists)
#   tier3_panel_conflicting.tsv      panel variants with conflicting ClinVar calls
#   germline_report.md               human-readable summary
#
# Interpretation rules of thumb (research-grade, NOT clinical):
#   * Tier 1 with >=2 stars in a panel gene -> take to a genetic counselor for a
#     clinical confirmatory test. Nothing here is actionable without that.
#   * Tier 2 are candidates only; most rare VUS are benign.
# =============================================================================
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
  library(data.table)
})

proj <- Sys.getenv("PROJ", getwd())   # run from the project folder, or set PROJ
out  <- file.path(proj, "results/01_germline")
rd   <- function(f) fread(file.path(out, f), sep = "\t", na.strings = c(".", "", "NA"), colClasses = "character")

exome <- rd("exome_clinvar_hits.tsv")
panel <- rd("panel_all_variants.tsv")
genes <- readLines(file.path(out, "panel_genes.txt"))
has_ann <- "ANN" %in% names(panel)

# ---------- helpers ----------
num <- function(x) suppressWarnings(as.numeric(x))
stars <- function(rev) fcase(
  grepl("practice_guideline", rev), 4L,
  grepl("reviewed_by_expert_panel", rev), 3L,
  grepl("multiple_submitters", rev) & grepl("no_conflicts", rev), 2L,
  grepl("single_submitter|conflicting", rev), 1L,
  default = 0L)
zyg <- function(gt) fcase(gt %in% c("1/1","1|1"), "hom",
                          gt %in% c("0/1","1/0","0|1","1|0"), "het",
                          grepl("^[0-9.]+[/|][0-9.]+$", gt), "other", default = NA_character_)
gene_from_geneinfo <- function(g) sub(":.*", "", sub("\\|.*", "", g))   # "TP53:7157|..." -> TP53

prep <- function(d) {
  d[, `:=`(zygosity = zyg(GT), stars = stars(CLNREVSTAT),
           gene = gene_from_geneinfo(GENEINFO),
           af = num(GNOMAD_AF), af_popmax = num(GNOMAD_AF_POPMAX),
           is_plp  = grepl("[Pp]athogenic", CLNSIG) & !grepl("Conflicting|Benign", CLNSIG),
           is_conf = grepl("Conflicting", CLNSIG))]
  if ("CALL_TYPE" %in% names(d))   # Sema4: GT is a bare "1", zygosity lives in INFO/CALL_TYPE
    d[is.na(zygosity), zygosity := fcase(CALL_TYPE == "het", "het", grepl("^hom", CALL_TYPE), "hom", default = NA_character_)]
  d[, alt_frac := vapply(strsplit(AD, ","), function(a) { a <- num(a); if (length(a) < 2 || sum(a) == 0) NA_real_ else a[2] / sum(a) }, 0)]
  if (has_ann) {   # SnpEff ANN: Allele|Annotation|Impact|Gene_Name|...  (first entry = most severe)
    a <- tstrsplit(sub(",.*", "", d$ANN), "|", fixed = TRUE)
    d[, `:=`(consequence = a[[2]], impact = a[[3]], ann_gene = a[[4]])]
    d[is.na(gene) | gene == "", gene := ann_gene]
  }
  d[]
}
exome <- prep(exome); panel <- prep(panel)

# gene symbol for panel variants that ClinVar doesn't know (no GENEINFO): look up from panel.bed
bed <- fread(file.path(out, "panel.bed"), col.names = c("chrom", "start", "end", "bed_gene"), colClasses = "character")
bed[, `:=`(start = as.integer(start), end = as.integer(end))]
panel[, posn := as.integer(POS)]
panel[is.na(gene) | gene == "", gene := bed[.SD, on = .(chrom = CHROM, start <= posn, end >= posn), x.bed_gene, mult = "first"]]
genes_covered <- unique(unlist(strsplit(bed$bed_gene, ",")))

# ---------- Tier 1: P/LP anywhere, >=1 star ----------
t1 <- exome[is_plp == TRUE & stars >= 1 & (is.na(af_popmax) | af_popmax < 0.05)]
t1[, in_panel := gene %in% genes]
setorder(t1, -in_panel, -stars)
# ClinVar "pathogenic" includes recessive carrier states and low-penetrance risk alleles;
# flag the ones that matter most: panel gene, or homozygous, or >=2 stars.
t1[, priority := fcase(in_panel & stars >= 2, "HIGH",
                       in_panel | zygosity == "hom" | stars >= 2, "MEDIUM", default = "LOW")]

# ---------- Tier 2: rare panel variants ----------
t2 <- panel[is_plp == FALSE & is_conf == FALSE & (is.na(af_popmax) | af_popmax < 0.005)]
t2 <- t2[!grepl("Benign", CLNSIG) | is.na(CLNSIG)]
if (has_ann) t2 <- t2[impact %in% c("HIGH", "MODERATE")]
t2[, gene_panel := gene %in% genes]
setorder(t2, -gene_panel, af_popmax, na.last = FALSE)

# ---------- Tier 3: conflicting in panel ----------
t3 <- panel[is_conf == TRUE]

# ---------- write ----------
keep <- intersect(c("CHROM","POS","REF","ALT","gene","zygosity","DP","AD","alt_frac","RED_FLAGS","CLNSIG","stars","CLNREVSTAT","CLNDN",
                    "af","af_popmax","GNOMAD_NHOMALT","consequence","impact","priority","in_panel","FILTER"), names(t1))
fwrite(t1[, ..keep], file.path(out, "tier1_pathogenic_exomewide.tsv"), sep = "\t")
keep2 <- intersect(c("CHROM","POS","REF","ALT","gene","zygosity","DP","AD","alt_frac","RED_FLAGS","CLNSIG","stars","af","af_popmax",
                     "consequence","impact","FILTER"), names(t2))
fwrite(t2[, ..keep2], file.path(out, "tier2_panel_rare.tsv"), sep = "\t")
fwrite(t3[, ..keep2], file.path(out, "tier3_panel_conflicting.tsv"), sep = "\t")

fmt_row <- function(d) if (nrow(d) == 0) "_none_" else d[, paste0("- **", gene, "** ", CHROM, ":", POS, " ", REF, ">", ALT,
  " (", zygosity, ", depth ", DP, sprintf(", alt fraction %.2f", alt_frac),
  if ("RED_FLAGS" %in% names(d)) ifelse(is.na(RED_FLAGS), "", paste0(", Sema4 flags: ", RED_FLAGS)) else "",
  ") — ", CLNSIG, " [", stars, "★]",
  ifelse(is.na(af_popmax), " gnomAD: n/a", sprintf(" gnomAD popmax AF %.2e", af_popmax)),
  if ("consequence" %in% names(d)) paste0(" — ", consequence) else "",
  ifelse(!is.na(CLNDN), paste0(" — ", substr(CLNDN, 1, 80)), ""))] |> paste(collapse = "\n")

rep <- c(
  "# Germline re-check — Sema4 normal (blood) exome vs current ClinVar/gnomAD",
  paste0("_generated ", Sys.Date(), "; hg19; research-grade only — confirm anything of interest with a clinical lab_"), "",
  "## Summary",
  sprintf("- ClinVar-annotated variants in exome: %d", nrow(exome)),
  sprintf("- Tier 1 (P/LP, ≥1 star, exome-wide): %d — HIGH: %d, MEDIUM: %d, LOW: %d",
          nrow(t1), sum(t1$priority == "HIGH"), sum(t1$priority == "MEDIUM"), sum(t1$priority == "LOW")),
  sprintf("- Tier 2 (rare, non-benign, in %d panel genes%s): %d", length(genes),
          if (has_ann) ", protein-altering" else ", ALL consequences — no ANN field; run VEP/OpenCRAVAT to narrow", nrow(t2)),
  sprintf("- Tier 3 (conflicting ClinVar, panel genes): %d", nrow(t3)),
  sprintf("- gnomAD frequencies available: %s", if (all(is.na(panel$af))) "NO (remote fetch failed; AF filters were skipped)" else "yes"), "",
  "## Tier 1 — HIGH / MEDIUM priority", fmt_row(t1[priority != "LOW"]), "",
  "## Tier 1 — LOW priority (carrier states, risk alleles, non-panel genes)", fmt_row(t1[priority == "LOW"]), "",
  "## Tier 2 — rare panel variants (candidates only)", fmt_row(t2[gene_panel == TRUE]), "",
  "## Tier 3 — conflicting classifications in panel genes", fmt_row(t3), "",
  "## Key negatives (report explicitly — a clean result is a finding)",
  paste0("- Panel genes with NO P/LP and NO rare protein-altering variant: ",
         paste(setdiff(intersect(genes, genes_covered), unique(c(t1[in_panel == TRUE]$gene, t2[gene_panel == TRUE]$gene, t3$gene))), collapse = ", ")),
  if (length(setdiff(genes, genes_covered))) paste0("- WARNING panel genes not found in refGene (not assessed): ", paste(setdiff(genes, genes_covered), collapse = ", ")) else "")
writeLines(rep, file.path(out, "germline_report.md"))
cat(rep, sep = "\n")
