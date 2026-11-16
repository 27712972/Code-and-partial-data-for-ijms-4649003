# Illustrative implementation; the original MR enrichment software is unconfirmed.
# Run from the repository root. Supply the retained MR list and/or module 09 output.
# Background: annotation-database genes, as specified in the latest Methods.
required <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors = FALSE, timeout = 600)
db <- org.Hs.eg.db::org.Hs.eg.db
inputs <- c(mr = "data/mr/prioritized_genes.txt",
            ko_down = "results/virtual_knockout/ko_down_genes.txt")
fdr_cutoff <- 0.05
top_n <- 5L
# These size limits follow the recovered KO enrichment script; MR settings unknown.
min_size <- 10L
max_size <- 500L

plot_terms <- function(d, path, color_limits, size_limits) {
  if (!nrow(d)) return(invisible(NULL))
  d$ratio <- vapply(strsplit(d$GeneRatio, "/", fixed = TRUE),
                    function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
  d$key <- paste(d$ontology, d$ID, sep = ":")
  labels <- setNames(d$Description, d$key)
  d$key <- factor(d$key, levels = rev(unique(d$key)))
  p <- ggplot(d, aes(ratio, key, size = Count, colour = -log10(pmax(p.adjust, 1e-300)))) +
    geom_point() + facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = labels) +
    scale_colour_gradient(low = "blue", high = "red", limits = color_limits) +
    scale_size_area(max_size = 8, limits = size_limits) +
    labs(x = "Gene ratio", y = NULL, colour = "-log10(BH FDR)", size = "Gene count") +
    theme_bw(base_size = 10)
  ggsave(path, p, width = 10, height = max(3, 0.35 * nrow(d) + 1.5), limitsize = FALSE)
}

run_set <- function(name, path) {
  genes <- unique(trimws(readLines(path, warn = FALSE)))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) stop("Empty gene list: ", path)
  out <- file.path("results/enrichment", name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(genes, file.path(out, "input_symbols.txt"))
  valid <- intersect(genes, AnnotationDbi::keys(db, keytype = "SYMBOL"))
  if (!length(valid)) stop("No valid human gene symbols in ", path)
  mapped <- AnnotationDbi::select(db, keys = valid, keytype = "SYMBOL", columns = "ENTREZID")
  mapping <- unique(merge(data.frame(SYMBOL = genes), mapped, all.x = TRUE, sort = FALSE))
  counts <- table(mapping$SYMBOL[!is.na(mapping$ENTREZID)])
  mapping$mapping_count <- as.integer(counts[mapping$SYMBOL])
  mapping$mapping_count[is.na(mapping$mapping_count)] <- 0L
  write.csv(mapping, file.path(out, "symbol_to_entrez_mapping.csv"), row.names = FALSE)
  ids <- unique(mapping$ENTREZID[!is.na(mapping$ENTREZID)])
  if (!length(ids)) stop("No Entrez identifiers mapped for ", name)
  writeLines(ids, file.path(out, "mapped_entrez_ids.txt"))
  summary_rows <- list()
  tops <- list()
  for (category in c("BP", "CC", "MF", "KEGG")) {
    # Omitting universe deliberately uses each annotation database's background.
    object <- if (category == "KEGG") {
      clusterProfiler::enrichKEGG(gene = ids, organism = "hsa", keyType = "ncbi-geneid",
        pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
        minGSSize = min_size, maxGSSize = max_size, use_internal_data = FALSE)
    } else {
      clusterProfiler::enrichGO(gene = ids, OrgDb = db, keyType = "ENTREZID", ont = category,
        pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1,
        minGSSize = min_size, maxGSSize = max_size, readable = TRUE)
    }
    prefix <- if (category == "KEGG") category else paste0("GO_", category)
    saveRDS(object, file.path(out, paste0(prefix, ".rds")))
    if (is.null(object)) {
      warning("No enrichment object returned for ", name, "/", category)
      next
    }
    all_terms <- as.data.frame(object@result)
    all_terms$ontology <- rep(category, nrow(all_terms))
    write.csv(all_terms, file.path(out, paste0(prefix, "_all_returned.csv")), row.names = FALSE)
    writeLines(as.character(object@universe), file.path(out, paste0(prefix, "_background.txt")))
    significant <- all_terms[!is.na(all_terms$p.adjust) & all_terms$p.adjust < fdr_cutoff, ]
    write.csv(significant, file.path(out, paste0(prefix, "_significant.csv")), row.names = FALSE)
    tops[[category]] <- head(significant[order(significant$p.adjust, significant$pvalue,
                                              seq_len(nrow(significant))), ], top_n)
    summary_rows[[category]] <- data.frame(analysis = category, input_symbols = length(genes),
      mapped_entrez = length(ids), background_genes = length(object@universe),
      returned_terms = nrow(all_terms), significant_terms = nrow(significant))
  }
  if (length(summary_rows)) write.csv(do.call(rbind, summary_rows),
                                    file.path(out, "analysis_summary.csv"), row.names = FALSE)
  displayed <- if (length(tops)) do.call(rbind, lapply(tops,
    function(d) d[, c("p.adjust", "Count"), drop = FALSE])) else data.frame()
  color_limits <- c(0, if (nrow(displayed))
    max(-log10(pmax(displayed$p.adjust, 1e-300))) else 1)
  size_limits <- c(0, if (nrow(displayed)) max(displayed$Count) else 1)
  for (category in intersect(c("BP", "KEGG"), names(tops))) {
    plot_terms(tops[[category]], file.path(out, paste0(category, "_top5.pdf")),
               color_limits, size_limits)
  }
  cc_mf <- tops[intersect(c("CC", "MF"), names(tops))]
  if (length(cc_mf)) plot_terms(do.call(rbind, cc_mf),
    file.path(out, "GO_CC_MF_top5_each.pdf"), color_limits, size_limits)
  writeLines(c("Background: annotation-database genes; no custom universe supplied.",
    "BH FDR < 0.05; top five significant terms per ontology by FDR, then P value.",
    "MR enrichment is a standardized demonstration; original software is unconfirmed.",
    "Live KEGG and installed annotation versions can change results.",
    capture.output(sessionInfo())), file.path(out, "run_notes.txt"))
}

available <- file.exists(inputs)
if (!any(available)) stop("No inputs found. Supply either: ", paste(inputs, collapse = " or "))
for (name in names(inputs)[!available]) message("Skipping missing input: ", inputs[[name]])
for (name in names(inputs)[available]) run_set(name, inputs[[name]])
