# Exploratory pooled-cell differential expression in mouse and human kidney.
# Dependencies: Seurat 5, ggplot2. Run from the package root after script 01.
# Input: processed RNA counts and log-normalized data, with group metadata.
# Output: all tested genes, selected genes, volcano plots and session information.
# These cell-level tests do not treat donors as biological replicates.
# min.pct = 0.25 is an illustrative choice from the available source script;
# it is not verified as the value that produced the final manuscript gene lists.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
out <- "results/differential_expression"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
min_pct <- 0.25
groups <- list(mouse = c("Con", "Mod"), human = c("Normal", "DKD"))

for (species in names(groups)) {
  input <- file.path("results/single_cell_processing", species,
                     paste0(species, "_kidney_processed.rds"))
  obj <- readRDS(input)
  stopifnot(inherits(obj, "Seurat"), "group" %in% names(obj[[]]))
  if (anyNA(obj$group) || !setequal(unique(obj$group), groups[[species]])) {
    stop("Expected both configured groups and no other group labels: ", species)
  }
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  normalized <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!nrow(normalized)) stop("Run RNA LogNormalize preprocessing first.")
  Idents(obj) <- "group"
  de <- FindMarkers(
    obj, ident.1 = groups[[species]][2], ident.2 = groups[[species]][1],
    assay = "RNA", slot = "data", test.use = "wilcox", min.pct = min_pct,
    logfc.threshold = 0, only.pos = FALSE, base = 2,
    fc.name = "avg_log2FC", pseudocount.use = 1, verbose = FALSE
  )
  # p_val_adj is Seurat's Bonferroni adjustment over all genes in the assay.
  result <- data.frame(gene = rownames(de), log2FC = de$avg_log2FC,
                       p_value = de$p_val, p_adj = de$p_val_adj)
  selected <- with(result, is.finite(log2FC) & !is.na(p_adj) &
                     p_adj < 0.01 & abs(log2FC) > 0.5)
  write.csv(result, file.path(out, paste0(species, "_all.csv")), row.names = FALSE)
  write.csv(result[selected, ], file.path(out, paste0(species, "_selected.csv")),
            row.names = FALSE)
  result$direction <- ifelse(selected,
    ifelse(result$log2FC > 0, "Higher in disease", "Lower in disease"),
    "Not selected")
  result$minus_log10_p <- -log10(pmax(result$p_adj, .Machine$double.xmin))
  plot <- ggplot(result, aes(log2FC, minus_log10_p, color = direction)) +
    geom_point(size = 0.7, alpha = 0.6, na.rm = TRUE) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(0.01), linetype = 2, color = "grey50") +
    scale_color_manual(values = c("Higher in disease" = "#C44E52",
      "Lower in disease" = "#4C72B0", "Not selected" = "grey75")) +
    labs(title = paste(species, "pooled-cell differential expression"),
         x = "Average log2 fold change: disease versus control",
         y = "-log10(Bonferroni-adjusted P)", color = NULL) + theme_bw()
  ggsave(file.path(out, paste0(species, "_volcano.pdf")), plot,
         width = 7, height = 5)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
