# Simple overview plots from analysis tables; not the final manuscript montage.
# Run from the project root. Each available module is plotted independently.
library(ggplot2)
out <- "results/summary_figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = width, height = height)
}
signed_colors <- c(down = "#3977B7", up = "#C64E52", other = "grey75")

# 1. Descriptive cell proportions, pooled within each group.
for (species in c("mouse", "human")) {
  path <- file.path("results/single_cell_processing", species, "cell_metadata.csv")
  if (!file.exists(path)) next
  d <- read.csv(path)
  if (all(d$celltype == "Unassigned")) {
    message("Skipping ", species, " composition: cell annotation is required.")
    next
  }
  d <- as.data.frame(table(group = d$group, celltype = d$celltype))
  d$proportion <- d$Freq / ave(d$Freq, d$group, FUN = sum)
  write.csv(d, file.path(out, paste0(species, "_cell_proportions.csv")), row.names = FALSE)
  p <- ggplot(d, aes(group, proportion, fill = celltype)) + geom_col() +
    labs(x = NULL, y = "Pooled cell proportion", fill = "Cell type") + theme_classic()
  save_plot(p, paste0(species, "_cell_proportions"), width = 9)
}

# 2. MR association counts, gene-cell matrix and RPS10 T-cell forest plot.
# The optional reported table must have gene, cell_type, b, se, pval and FDR.
# It may be a curated display subset: never recompute its FDR locally.
mr_file <- "data/mr/reported_results.csv"
if (!file.exists(mr_file)) mr_file <- "results/mr/primary_results.csv"
if (file.exists(mr_file)) {
  all_mr <- read.csv(mr_file, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "cell_type", "b", "se", "pval", "FDR") %in% names(all_mr)),
            !anyDuplicated(all_mr[c("gene", "cell_type")]))
  d <- subset(all_mr, is.finite(FDR) & FDR >= 0 & FDR < 0.05 & is.finite(b))
  if (nrow(d)) {
    d$direction <- ifelse(d$b < 0, "down", "up")
    p <- ggplot(d, aes(cell_type, fill = direction)) + geom_bar() + coord_flip() +
      scale_fill_manual(values = signed_colors, labels = c(down = "OR < 1", up = "OR >= 1")) +
      labs(x = NULL, y = "Selected gene-cell associations", fill = NULL) + theme_classic()
    save_plot(p, "MR_association_counts", height = 6)
    d$significance <- -log10(pmax(d$FDR, .Machine$double.xmin))
    p <- ggplot(d, aes(gene, cell_type, color = b, size = significance)) + geom_point() +
      scale_color_gradient2(low = "#3977B7", mid = "white", high = "#C64E52") +
      scale_size_area(max_size = 5) + labs(x = NULL, y = NULL, color = "log(OR)",
      size = "-log10(FDR)") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
    save_plot(p, "MR_gene_cell_matrix", width = 15, height = 7)
  }
  t_types <- c("CD4 Effector memory/TEMRA", "CD4 Naive/Central memory T cell",
               "CD8 Effector memory", "CD8 Naive/Central memory T cell")
  d <- subset(all_mr, gene == "RPS10" & cell_type %in% t_types &
                is.finite(b) & is.finite(se) & se > 0)
  if (nrow(d)) {
    p <- ggplot(d, aes(cell_type, exp(b))) + geom_hline(yintercept = 1, linetype = 2) +
      geom_errorbar(aes(ymin = exp(b - 1.96 * se), ymax = exp(b + 1.96 * se)),
                    width = 0.15) + geom_point(color = "#C64E52", size = 3) +
      coord_flip() + scale_y_log10() + labs(x = NULL, y = "Odds ratio (95% CI)") +
      theme_classic()
    save_plot(p, "RPS10_Tcell_MR_forest")
  }
}

# 3. Ribosomal-exclusion knockout predictions and shared-gene rank agreement.
ko_file <- "results/virtual_knockout/DKD/remove_other_ribosomal_KO.csv"
if (file.exists(ko_file)) {
  d <- read.csv(ko_file)
  d <- subset(d, !is_target & is.finite(direction_score_depth_adjusted) &
                is.finite(network_FDR_BH))
  d$display <- ifelse(d$selected, d$predicted_direction_depth_adjusted, "other")
  d$display[!d$display %in% c("up", "down")] <- "other"
  p <- ggplot(d, aes(direction_score_depth_adjusted,
                     -log10(pmax(network_FDR_BH, .Machine$double.xmin)), color = display)) +
    geom_point(size = 1.3) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
    scale_color_manual(values = signed_colors) + labs(x = "Depth-adjusted direction score",
      y = "-log10(model FDR)", color = NULL) + theme_classic()
  save_plot(p, "RPS10_KO_direction")
}
rank_file <- "results/virtual_knockout/DKD/shared_gene_rank_comparison.csv"
if (file.exists(rank_file)) {
  d <- read.csv(rank_file)
  d$display <- "other"
  top <- head(order(-d$exclusion_distance, d$gene), 20)
  d$display[top] <- d$predicted_direction[top]
  d$display[!d$display %in% c("up", "down")] <- "other"
  p <- ggplot(d, aes(primary_rank_fraction, exclusion_rank_fraction, color = display)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 1.3) + scale_color_manual(values = signed_colors) + coord_equal() +
    labs(x = "Primary perturbation-distance rank fraction",
         y = "Ribosomal-exclusion rank fraction", color = NULL) + theme_classic()
  save_plot(p, "RPS10_KO_rank_agreement")
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
