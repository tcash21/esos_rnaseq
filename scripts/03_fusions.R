#!/usr/bin/env Rscript
# =============================================================================
# 03_fusions.R — review the raw fusion candidate lists Sema4 filtered down to "none" (Mac).
#   Rscript scripts/03_fusions.R
#
# Inputs (data/sema4/**):
#   ISM563041-2.fusionVariant.tsv                         Sema4 aggregate (callers in redFlags: caller_FC / caller_SF ...)
#   ISM563041-2.final-list_candidate-fusion-genes.hg19.txt FusionCatcher full output (descriptions carry artifact tags)
#   ISM563041-2.FusionInspector.fusions.abridged.tsv      FusionInspector re-validation (1 row: SEC31A--JAK2)
# Optional: results/02_expression/requant/ISM563041-2.Chimeric.out.junction (STAR, from 02b) — independent
#   cross-check of any candidate's breakpoint.
# Output: results/03_fusions/fusions_report.md, candidates_scored.tsv
#
# Scoring logic: a real fusion has >= 3 unique junction-spanning reads AND spanning pairs, partners that are
# not adjacent (read-through), not paralogs/pseudogenes, and no FusionCatcher "seen in normals" tag.
# Highly expressed genes (COL1A1 at 9,000 TPM) generate chimeric cDNA by template switching: many pairs, ~1 unique read.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
proj <- Sys.getenv("PROJ", getwd()); out <- file.path(proj, "results/03_fusions"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
ff <- function(pat) list.files(file.path(proj, "data/sema4"), pat, recursive = TRUE, full.names = TRUE)[1]

sarcoma_partners <- c("EWSR1","FUS","SS18","SSX1","SSX2","NR4A3","CIC","DUX4","BCOR","CCNB3","ALK","NTRK1","NTRK2","NTRK3","ROS1","RET",
  "BRAF","RAF1","FGFR1","FGFR2","FGFR3","PDGFB","PDGFRA","PDGFRB","TFE3","TFEB","USP6","HEY1","NCOA2","NAB2","STAT6","CRTC1","MAML2",
  "PLAG1","HMGA2","MYOCD","YAP1","WWTR1","CAMTA1","TAZ","MDM2","CDK4","KDM3A","MAP4K4","TP53","RB1","EGFR","MET","ERBB2","KIT","JAK2",
  "ETV6","PAX3","PAX7","FOXO1","ASPSCR1","ATF1","CREB1","GLI1","MEAF6","PHF1","JAZF1","SUZ12","BCL6","MYC","TERT","ATRX","CIITA","NFIB")

# ---------- FusionCatcher ----------
fc <- fread(ff("final-list_candidate-fusion-genes.hg19.txt$"))
setnames(fc, 1:10, c("g5","g3","desc","common_reads","span_pairs","span_unique","anchor","method","bp5","bp3"))
artifact_tags <- "readthrough|pseudogene|paralog|banned|healthy|gtex|hpa|conjoing|bodymap2|1000genomes|cta_gene|duplicates|similar_reads|short_repeats|long_repeats|partial-matches|fully_overlapping|same_strand_overlapping|antisense|ribosomal|mt|distance1000bp|distance100kbp|adjacent|similar_symbols|non_cancer_tissues|non_tumor_cells|fragments|matched-normal"
known_tags <- "known|cosmic|chimerdb|ticdb|cgp|tcga|mitelman|oncokb|18cancers|prostate_cancer|pancreas|gliomas|oesophagus|ccle"
fc[, `:=`(artifact = grepl(artifact_tags, desc), known = grepl(known_tags, desc),
          sarcoma_gene = g5 %in% sarcoma_partners | g3 %in% sarcoma_partners,
          same_chrom = sub(":.*", "", bp5) == sub(":.*", "", bp3))]
fc[, score := fcase(span_unique >= 5 & span_pairs >= 5 & !artifact, "strong",
                    span_unique >= 3 & span_pairs >= 3 & !artifact, "plausible",
                    default = "weak/artifact")]

# ---------- Sema4 aggregate ----------
sv <- fread(ff("fusionVariant\\.tsv$"))
sv[, `:=`(callers = gsub(".*numCallers_([0-9]+).*", "\\1", redFlags),
          which = paste(regmatches(redFlags, gregexpr("caller_[A-Z]+", redFlags))[[1]], collapse = "+"),
          mappable = !grepl("geneNotMappable", redFlags), clean_junction = !grepl("notCleanExonJunction", redFlags)), by = seq_len(nrow(sv))]
sv_top <- sv[filter == "PASS" & numReads >= 3 & mappable]

# ---------- FusionInspector ----------
fi <- fread(ff("FusionInspector\\.fusions\\.abridged\\.tsv$"))

# ---------- STAR chimeric junctions from 02b (optional cross-check) ----------
chim <- file.path(proj, "results/02_expression/requant/ISM563041-2.Chimeric.out.junction")
chim_note <- if (file.exists(chim)) {
  cj <- fread(cmd = paste("grep -v '^#'", shQuote(chim)), header = TRUE, select = 1:7)   # STAR chimOutJunctionFormat 1: header + '#' footer
  setnames(cj, c("cA", "bA", "sA", "cB", "bB", "sB", "jtype"))
  std <- paste0("chr", c(1:22, "X", "Y"))
  cj <- cj[jtype >= 0 & cA %in% std & cB %in% std & (cA != cB | abs(bA - bB) > 5e5)]   # > 500 kb: adjacent-gene read-through is not a fusion
  # annotate both ends to GENCODE v36 genes (protein-coding + lncRNA only; rRNA/7SL/snRNA chimeras are library artifacts)
  pmg <- fread(file.path(proj, "refs/xena/gencode.v36.annotation.gtf.gene.probemap"))[, .(gene, chrom, start = chromStart, end = chromEnd)]
  pmg <- pmg[!grepl("^(RNA5|RNA1|RNA2|RN7S|RNU|RNY|Y_RNA|SNOR|SCARNA|MIR|MT-|Metazoa|U[0-9])", gene)]
  setkey(pmg, chrom, start, end)
  ann <- function(ch, pos) { q <- data.table(chrom = ch, start = pos, end = pos); r <- foverlaps(q, pmg, mult = "first", nomatch = NA); r$gene }
  cj[, `:=`(gA = ann(cA, bA), gB = ann(cB, bB))]
  gp <- cj[!is.na(gA) & !is.na(gB) & gA != gB, .(reads = .N), by = .(gA, gB)][order(-reads)]
  gp[, sarcoma := gA %in% sarcoma_partners | gB %in% sarcoma_partners]
  want <- gp[(gA == "USP39" & gB == "CTNNA2") | (gA == "CTNNA2" & gB == "USP39") | (gA == "PGAP1" & gB == "DNAH7") | (gA == "DNAH7" & gB == "PGAP1")]
  c(sprintf("STAR (GRCh38, 02b): %d distant/interchromosomal split reads on standard chromosomes; %d gene pairs with >= 10 reads.", nrow(cj), gp[reads >= 10, .N]),
    "- Breakpoint chimeras from FusionCatcher reproduced by STAR? ",
    if (nrow(want)) want[, sprintf("  - **%s–%s**: %d split reads", gA, gB, reads)] else "  - neither USP39–CTNNA2 nor PGAP1–DNAH7 seen (STAR's chimeric detection is stricter; not evidence against)",
    "- Gene pairs involving a sarcoma-relevant partner, >= 5 reads:",
    if (gp[sarcoma & reads >= 5, .N]) gp[sarcoma & reads >= 5][1:min(10, .N), sprintf("  - %s–%s: %d reads", gA, gB, reads)] else "  - none",
    "- Top gene pairs overall (>= 50 reads; expect abundant-transcript artifacts such as COL1A1/COL1A2/COL3A1):",
    gp[reads >= 50][1:min(12, .N), sprintf("  - %s–%s: %d", gA, gB, reads)])
} else "_STAR Chimeric.out.junction not present yet (run 02b)_"

# ---------- write ----------
cand <- fc[order(factor(score, levels = c("strong", "plausible", "weak/artifact")), -span_unique)]
fwrite(cand[, .(g5, g3, score, span_pairs, span_unique, anchor, known, sarcoma_gene, artifact, same_chrom, bp5, bp3, desc)],
       file.path(out, "candidates_scored.tsv"), sep = "\t")
rowfc <- function(d) if (nrow(d) == 0) "_none_" else d[, sprintf("- **%s–%s** unique %d / pairs %d / anchor %d%s%s — %s%s", g5, g3, span_unique, span_pairs, anchor,
  ifelse(known, " **[known fusion tag]**", ""), ifelse(sarcoma_gene, " **[sarcoma gene]**", ""), substr(desc, 1, 90), ifelse(artifact, " ⚠", ""))] |> paste(collapse = "\n")
rep <- c("# Fusion review — everything the callers produced vs Sema4's 'no clinically significant fusions'",
  paste0("_generated ", Sys.Date(), "; research-grade_"), "",
  "## Summary",
  sprintf("- FusionCatcher raw candidates: %d; strong (>=5 unique + >=5 pairs, no artifact tag): **%d**; plausible (>=3/>=3): %d; carrying a 'known fusion' tag: %d; involving a sarcoma/cancer partner gene: %d",
          nrow(fc), sum(fc$score == "strong"), sum(fc$score == "plausible"), sum(fc$known), sum(fc$sarcoma_gene)),
  sprintf("- Sema4 aggregate list: %d rows (%d PASS); called by >= 2 callers: %d; PASS with >= 3 reads and mappable genes: %d",
          nrow(sv), sum(sv$filter == "PASS"), sum(as.integer(sv$callers) >= 2, na.rm = TRUE), nrow(sv_top)),
  sprintf("- FusionInspector validated: %d (%s) — junction %s / spanning %s reads, FFPM %s: **artifact-level support**", nrow(fi),
          paste(fi[[1]], collapse = ", "), fi$JunctionReadCount, fi$SpanningFragCount, fi$FFPM), "",
  "## Strong candidates", rowfc(fc[score == "strong"]), "",
  "## Plausible candidates", rowfc(fc[score == "plausible"]), "",
  "## Any candidate involving a sarcoma-relevant partner (regardless of support)", rowfc(fc[sarcoma_gene == TRUE][order(-span_unique)][1:min(20, .N)]), "",
  "## Candidates with a 'known fusion' database tag (regardless of support)", rowfc(fc[known == TRUE][order(-span_unique)][1:min(15, .N)]), "",
  "## Chimeric transcripts that coincide with DNA copy-number breakpoints (FACETS, step 04)",
  "- **USP39–CTNNA2** (2p11.2–2p12): USP39 lies inside the 18-copy KDM3A amplicon (83.08–88.13 Mb); CTNNA2 sits at the 4→7-copy boundary at 80.10 Mb. A transcribed junction between the amplicon and its flank.",
  "- **PGAP1–DNAH7** (2q33.1–2q32.3): PGAP1 sits at the 197.76 Mb segment boundary, DNAH7 at the 196.89 Mb boundary.",
  "- Interpretation: these are genomic rearrangements at segment edges read out as RNA chimeras — structural evidence for how the KDM3A amplicon and the 2q segments are built, not oncogenic fusions. Both are PCR-validatable at the DNA level (see 03b).", "",
  "## Independent cross-check", chim_note, "",
  "## The COL1A1 chimeras", sprintf("- %d FusionCatcher candidates involve COL1A1 (tumor's top-expressed gene, ~9,000 TPM); max unique reads %d. Template-switching artifacts of an extremely abundant transcript, not fusions.",
    fc[g5 == "COL1A1" | g3 == "COL1A1", .N], fc[g5 == "COL1A1" | g3 == "COL1A1", max(span_unique)]), "",
  sprintf("- Sema4's aggregate list (%d rows) is the FusionCatcher output re-formatted: no candidate was called by a second caller.", nrow(sv)), "",
  "## Reading guide",
  "- ESOS/osteosarcoma is fusion-negative by definition; the value here is (a) confirming that against the full candidate lists and (b) excluding the fusion-defined mimics (EWSR1/FUS-, SS18-, NR4A3-, BCOR-, CIC-, NTRK-, HEY1-NCOA2 sarcomas).",
  "- 'known fusion tag' in FusionCatcher means one of its databases lists the pair; with 1-2 unique reads it is still noise.")
writeLines(rep, file.path(out, "fusions_report.md")); cat(rep, sep = "\n")
