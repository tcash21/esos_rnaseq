#!/usr/bin/env Rscript
# =============================================================================
# 05_somatic.R — triage Sema4's somatic calls after VEP (Mac).   Rscript scripts/05_somatic.R
#
# Inputs : results/05_somatic/somatic_annotated.tsv      (05_somatic.sh, EC2)
#          results/04_copynumber/facets_cval150_{fit,segments}.tsv   (purity + local copy number)
# Outputs: results/05_somatic/
#   somatic_clean.tsv      calls passing the filter policy below, with CCF and multiplicity
#   somatic_report.md
#
# Filter policy (see NOTEBOOK §4 item 6): Sema4's PASS set is dominated by calls Mutect2 had
# filtered and Sema4's pipeline rescued ("mutectFiltOverride"), nearly all below 5% AF.
# We keep: FILTER==PASS, >= 8 alt reads, tumor AF >= 0.05 (SNV) / >= 0.10 (indel), normal AF <= 0.02,
# and drop techFail/badDp. Sema4's lowAfT / readsInN flags are replaced by these explicit thresholds
# (readsInN fires on a single normal read). "mutectFiltOverride" is NOT excluded: Mutect2's
# clustered-events filter removes true kataegis clusters (MAP4K4 here). Calls within 10 bp of each
# other are flagged as one event (indel realignment ambiguity / read-cluster noise).
# =============================================================================
suppressPackageStartupMessages(library(data.table))
proj <- Sys.getenv("PROJ", getwd())
out  <- file.path(proj, "results/05_somatic"); cn <- file.path(proj, "results/04_copynumber")
num  <- function(x) suppressWarnings(as.numeric(x))

d <- fread(file.path(out, "somatic_annotated.tsv"), na.strings = c(".", ""), colClasses = "character")
d[, `:=`(t_af = num(T_AF), t_dp = num(T_DP), n_af = num(N_AF), n_dp = num(N_DP),
         t_alt = num(sub(".*,", "", T_AD)), CHROM = as.character(CHROM), POS = as.integer(POS))]
d[, hgvsp := sub("%3D", "=", sub("^[^:]*:", "", HGVSp), fixed = TRUE)]
d[, hgvsc := sub("^[^:]*:", "", HGVSc)]
d[, bad_flag := grepl("techFail|badDp", RED_FLAGS)]
d[, is_snv := nchar(REF) == 1 & nchar(ALT) == 1]
d[, keep := FILTER == "PASS" & !bad_flag & t_alt >= 8 & t_af >= ifelse(is_snv, 0.05, 0.10) & (is.na(n_af) | n_af <= 0.02)]
setorder(d, CHROM, POS)
d[, near := c(FALSE, diff(POS) <= 10 & CHROM[-1] == CHROM[-.N]) | c(diff(POS) <= 10 & CHROM[-1] == CHROM[-.N], FALSE)]

fit  <- fread(file.path(cn, "facets_cval150_fit.tsv")); purity <- fit$purity
segs <- fread(file.path(cn, "facets_cval150_segments.tsv"))[, chrom := as.character(chrom)]
d <- segs[d, on = .(chrom = CHROM, start <= POS, end >= POS), mult = "first",
          .(CHROM = i.CHROM, POS = i.POS, REF, ALT, FILTER, RED_FLAGS, SYMBOL, Consequence, IMPACT, hgvsc, hgvsp, EXON,
            SIFT, PolyPhen, Existing_variation, MAX_AF, t_af, t_dp, t_alt, n_af, n_dp, keep, is_snv, near, tcn = x.tcn.em, lcn = x.lcn.em)]
# multiplicity (mutated copies per tumor cell) and CCF (Dentro 2017):  AF = p*m / (p*tcn + 2(1-p)) * ccf
d[, mult := t_af * (purity * tcn + 2 * (1 - purity)) / purity]
d[, `:=`(ccf = pmin(1, mult / pmax(1, round(mult))), mult_round = pmax(1, round(mult)))]
d[, timing := fcase(is.na(tcn), NA_character_,
                    tcn >= 3 & mult_round >= 2, "early (before gain/WGD)",
                    tcn >= 3 & mult_round == 1, "late (after gain/WGD) or subclonal",
                    default = "n/a (<=2 copies)")]

clean <- d[keep == TRUE][order(-t_af)]
fwrite(clean, file.path(out, "somatic_clean.tsv"), sep = "\t")

# kataegis-like clusters: >= 3 kept SNVs within 2 kb
cs <- clean[is_snv == TRUE][order(CHROM, POS)]
cs[, grp := cumsum(c(TRUE, diff(POS) > 2000)), by = CHROM]
cl <- cs[, .(n = .N, span = max(POS) - min(POS), genes = paste(unique(SYMBOL), collapse = ","), af = sprintf("%.2f", mean(t_af)),
      subs = paste(sort(unique(paste0(REF, ">", ALT))), collapse = ",")), by = .(CHROM, grp)][n >= 3]

spec <- clean[is_snv == TRUE, .(sub = fifelse(REF %in% c("C", "T"), paste0(REF, ">", ALT),
        paste0(chartr("ACGT", "TGCA", REF), ">", chartr("ACGT", "TGCA", ALT))))][, .N, by = sub][order(sub)]

tmb <- clean[IMPACT %in% c("HIGH", "MODERATE"), .N] / 35   # ~35 Mb capture, nonsynonymous per Mb (Sema4 definition)
cancer <- c("TP53","RB1","ATRX","BRCA1","BRCA2","PALB2","RAD51C","RAD51D","BRIP1","ATM","CHEK2","PTEN","NF1","CDKN2A","MDM2","CDK4",
            "KDM3A","MAP4K4","DLG2","PIK3CA","KRAS","NRAS","HRAS","IDH1","IDH2","SETD2","ARID1A","KMT2C","KMT2D","CREBBP","EP300",
            "H3F3A","H3F3B","PRKDC","DICER1","MYC","CCNE1","AKT1","LRP1B","FANCA","FANCD2","FANCM","POLE","SMARCA4","SMARCB1",
            "TERT","APC","CTNNB1","FGFR1","STAG2","KEAP1","NFE2L2","ARID1B","BCOR","ZFHX3")
f2 <- function(x) sprintf("%.2f", x)
row <- function(x) x[, sprintf("- **%s** %s %s (%s) — AF %.2f, %d/%d reads, normal %d/%d; CN %s/%s, mult %s, CCF %s — %s%s%s",
        SYMBOL, ifelse(is.na(hgvsp), hgvsc, hgvsp), Consequence, IMPACT, t_af, t_alt, t_dp, round(n_af * n_dp), n_dp,
        ifelse(is.na(tcn), "?", tcn), ifelse(is.na(lcn), "?", lcn), ifelse(is.na(mult), "?", f2(mult)), ifelse(is.na(ccf), "?", f2(ccf)),
        ifelse(is.na(timing), "", timing), ifelse(is.na(Existing_variation), "", paste0(" [", Existing_variation, "]")),
        ifelse(near, " ⚠ adjacent call within 10 bp: likely one event", ""))] |> paste(collapse = "\n")

rep <- c("# Somatic re-annotation — Sema4 calls vs VEP 116 + FACETS", paste0("_generated ", Sys.Date(), "; research-grade_"), "",
  "## Filter policy result",
  sprintf("- Sema4 PASS calls: %d; of these %d carry mutectFiltOverride (Mutect2 had filtered them), %d lowAfT, %d readsInN",
          d[FILTER == "PASS", .N], d[FILTER == "PASS" & grepl("mutectFiltOverride", RED_FLAGS), .N],
          d[FILTER == "PASS" & grepl("lowAfT", RED_FLAGS), .N], d[FILTER == "PASS" & grepl("readsInN", RED_FLAGS), .N]),
  sprintf("- Kept after policy (PASS, >=8 alt reads, AF>=0.05 SNV / >=0.10 indel, normal AF<=0.02): **%d** (%d SNV, %d indel; %d flagged as adjacent-call pairs)",
          nrow(clean), clean[is_snv == TRUE, .N], clean[is_snv == FALSE, .N], clean[near == TRUE, .N]),
  sprintf("- Nonsynonymous TMB from the kept set: %.2f /Mb (Sema4 reported 3.77; definition of their denominator/numerator unknown)", tmb),
  sprintf("- Purity used for CCF/multiplicity: %.2f (FACETS)", purity), "",
  "## Kept calls in cancer-relevant genes", row(clean[SYMBOL %in% cancer & IMPACT %in% c("HIGH", "MODERATE")]), "",
  "## Kataegis-like clusters (>=3 kept SNVs within 2 kb)",
  if (nrow(cl)) cl[, sprintf("- chr%s: %d SNVs over %d bp in %s, mean AF %s, substitutions %s", CHROM, n, span, genes, af, subs)] else "_none_", "",
  "## Substitution spectrum of kept SNVs (pyrimidine context)",
  paste(spec[, sprintf("%s: %d", sub, N)], collapse = "; "), "",
  "## Clonality of kept calls",
  sprintf("- CCF >= 0.8 (clonal): %d; CCF < 0.8: %d; multiplicity >= 2 in gained regions (pre-WGD): %d",
          clean[ccf >= 0.8, .N], clean[ccf < 0.8, .N], clean[timing == "early (before gain/WGD)", .N]), "",
  "## All kept HIGH-impact calls", row(clean[IMPACT == "HIGH"]), "",
  "## Discarded but worth knowing (PASS in Sema4, failed policy) — cancer genes only",
  row(d[keep == FALSE & FILTER == "PASS" & SYMBOL %in% cancer & IMPACT %in% c("HIGH", "MODERATE")][order(-t_af)][1:min(15, .N)]),
  "  _(nearly all <5% AF with Mutect2-override flags: read-cluster noise, not mutations)_")
writeLines(rep, file.path(out, "somatic_report.md")); cat(rep, sep = "\n")
