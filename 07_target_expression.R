# Descriptive cell-level expression of RPL8/RPS10 and their mouse orthologues.
# Dependencies: Seurat 5, ggplot2. Run from the package root after script 01.
# Input: annotated processed objects; RNA data must be log-normalized.
# Output: cell-level values, unadjusted Wilcoxon P values and expression plots.
# Cells are correlated within donors; these are exploratory cell-level tests,
# not donor-level evidence. Script 08 performs the donor-level analysis.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
human_types <- c("B cells", "T cells", "Dendritic cells", "Monocytes",
                 "Macrophages", "Mast cells", "Neutrophils", "pDC")
# Example spelling: replace only with the verified broad mouse immune label.
mouse_immune_label <- "Immune cells"
config <- list(
  mouse = list(genes = c("Rpl8", "Rps10"), groups = c("Con", "Mod"),
               types = mouse_immune_label),
  human = list(genes = c("RPL8", "RPS10"), groups = c("Normal", "DKD"),
               types = human_types)
)
out <- "results/target_expression"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (species in names(config)) {
  cfg <- config[[species]]
  obj <- readRDS(file.path("results/single_cell_processing", species,
                          paste0(species, "_kidney_processed.rds")))
  stopifnot(inherits(obj, "Seurat"))
  md <- obj[[]]
  stopifnot(all(c("donor", "group", "celltype") %in% names(md)))
  md[c("donor", "group", "celltype")] <-
    lapply(md[c("donor", "group", "celltype")], as.character)
  if (anyNA(md[, c("donor", "group", "celltype")])) stop("Missing metadata.")
  if (!setequal(unique(md$group), cfg$groups)) stop("Unexpected group labels.")
  if (!any(md$celltype %in% cfg$types)) {
    stop("No configured immune labels found. Complete and check annotation: ", species)
  }
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  x <- GetAssayData(obj, assay = "RNA", layer = "data")
  if (!all(cfg$genes %in% rownames(x))) stop("Target genes missing: ", species)
  cells <- rownames(md)[md$celltype %in% cfg$types]
  md <- md[cells, , drop = FALSE]
  values <- as.matrix(x[cfg$genes, cells, drop = FALSE])
  long <- do.call(rbind, lapply(cfg$genes, function(gene) data.frame(
    cell_id = cells, donor = md$donor, group = md$group,
    celltype = md$celltype, gene = gene, expression = as.numeric(values[gene, ])
  )))
  long$group <- factor(long$group, levels = cfg$groups)
  long$celltype <- factor(long$celltype, levels = cfg$types)
  tests <- expand.grid(gene = cfg$genes, celltype = cfg$types,
                       stringsAsFactors = FALSE)
  tests$n_control <- tests$n_disease <- 0L
  tests$p_value <- NA_real_
  tests$status <- "Insufficient cells"
  for (i in seq_len(nrow(tests))) {
    d <- long[long$gene == tests$gene[i] & long$celltype == tests$celltype[i], ]
    n <- table(d$group)
    tests$n_control[i] <- n[1]
    tests$n_disease[i] <- n[2]
    # Three cells per group is an example guard, not a biological replicate rule.
    if (any(n < 3L)) next
    test <- pairwise.wilcox.test(d$expression, d$group, paired = FALSE,
      alternative = "two.sided", p.adjust.method = "none", exact = FALSE,
      correct = TRUE)
    p <- as.numeric(test$p.value[1, 1])
    tests$p_value[i] <- if (is.finite(p)) p else NA_real_
    tests$status[i] <- if (is.finite(p)) "Tested" else "Undefined statistic"
  }
  write.csv(long, file.path(out, paste0(species, "_cell_expression.csv")),
            row.names = FALSE)
  write.csv(tests, file.path(out, paste0(species, "_wilcoxon.csv")),
            row.names = FALSE)
  tests$label <- ifelse(is.na(tests$p_value), "Not tested",
                        paste0("P = ", format.pval(tests$p_value, digits = 2)))
  tests$celltype <- factor(tests$celltype, levels = cfg$types)
  p <- ggplot(long, aes(group, expression, fill = group)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
    geom_text(data = tests, aes(x = 1.5, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.3, size = 2.6) +
    facet_grid(gene ~ celltype, drop = FALSE) +
    scale_fill_manual(values = c("#4C72B0", "#C44E52")) +
    labs(x = NULL, y = "RNA log-normalized expression", fill = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        legend.position = "none")
  ggsave(file.path(out, paste0(species, "_expression.pdf")), p,
         width = max(5, 2 * length(cfg$types)), height = 6)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
