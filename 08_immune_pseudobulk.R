# Donor-level human immune pseudobulk: RPL8/RPS10, eight types, 16 hypotheses.
# Dependencies: Seurat 5, Matrix, edgeR with modern QL support, statmod, ggplot2.
# Input: annotated human object from script 01; donor identifies biological donors.
# Output: raw sums, coverage, model status, target estimates, LODO and overview plots.
# Counts are summed without log transformation or division by cell number.
# Main eligibility: >=20 cells per profile and >=3 donors in each group.
# The 100-gene and two-residual-df safeguards follow the available source code.
# LODO permits two donors in the omitted donor's group: sensitivity, not replication.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})
if (!requireNamespace("statmod", quietly = TRUE)) stop("Install statmod first.")
if (!"legacy" %in% names(formals(getFromNamespace("glmQLFit.default", "edgeR")))) {
  stop("An edgeR version supporting modern QL (legacy = FALSE) is required.")
}
immune_types <- c("B cells", "T cells", "Dendritic cells", "Monocytes",
                 "Macrophages", "Mast cells", "Neutrophils", "pDC")
targets <- c("RPL8", "RPS10")
out <- "results/immune_pseudobulk"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS("results/single_cell_processing/human/human_kidney_processed.rds")
stopifnot(inherits(obj, "Seurat"))
if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
md <- obj[[]][colnames(counts), , drop = FALSE]
stopifnot(all(c("donor", "group", "celltype") %in% names(md)))
md[c("donor", "group", "celltype")] <-
  lapply(md[c("donor", "group", "celltype")], as.character)
if (anyNA(md[, c("donor", "group", "celltype")])) stop("Missing metadata.")
if (!setequal(unique(md$group), c("Normal", "DKD"))) stop("Unexpected groups.")
if (!any(md$celltype %in% immune_types)) stop("Complete the human cell annotation.")
v <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
if (any(!is.finite(v) | v < 0 | abs(v - round(v)) > 1e-7)) {
  stop("Use original, nonnegative integer RNA counts.")
}
donors <- unique(md[, c("donor", "group")])
if (anyDuplicated(donors$donor)) stop("A biological donor cannot have two groups.")
if (any(!nzchar(trimws(donors$donor)))) stop("Blank donor identifiers.")
coverage <- expand.grid(donor = donors$donor, celltype = immune_types,
                        stringsAsFactors = FALSE)
coverage$group <- donors$group[match(coverage$donor, donors$donor)]
coverage$n_cells <- as.vector(table(factor(md$donor, levels = donors$donor),
                                   factor(md$celltype, levels = immune_types)))
coverage$eligible <- coverage$n_cells >= 20L
write.csv(coverage, file.path(out, "donor_celltype_coverage.csv"), row.names = FALSE)

fit_ql <- function(x, group, retained = NULL) {
  group <- factor(group, levels = c("Normal", "DKD"))
  design <- model.matrix(~ group)
  if (qr(design)$rank < ncol(design) || nrow(design) - ncol(design) < 2L) {
    stop("The group design has insufficient residual degrees of freedom.")
  }
  y <- DGEList(counts = as.matrix(x))
  keep <- if (is.null(retained)) filterByExpr(y, design = design) else
    rownames(y) %in% retained & rowSums(y$counts) > 0
  if (sum(keep) < 100L) stop("Fewer than 100 expressed genes for model fitting.")
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (any(y$samples$lib.size <= 0)) stop("Empty library after gene filtering.")
  y <- calcNormFactors(y, method = "TMM")
  fit <- glmQLFit(y, design = design, robust = TRUE, legacy = FALSE)
  test <- glmQLFTest(fit, coef = "groupDKD")
  list(table = test$table, genes = rownames(y), normalized = y)
}

result <- expand.grid(gene = targets, celltype = immune_types, stringsAsFactors = FALSE)
result$log2FC <- result$lnFC <- result$p_value <- NA_real_
result$status <- "Insufficient donor coverage"
audit <- data.frame(celltype = immune_types, n_control = 0L, n_disease = 0L,
                     status = "Insufficient donor coverage")
raw_sums <- models <- lodo <- list()
for (ct in immune_types) {
  cells <- which(md$celltype == ct)
  membership <- sparseMatrix(i = seq_along(cells),
    j = match(md$donor[cells], donors$donor), x = rep(1, length(cells)),
    dims = c(length(cells), nrow(donors)))
  pb <- counts[, cells, drop = FALSE] %*% membership
  colnames(pb) <- donors$donor
  raw_sums[[ct]] <- pb
  eligible <- coverage[coverage$celltype == ct & coverage$eligible, ]
  n <- table(factor(eligible$group, levels = c("Normal", "DKD")))
  a <- match(ct, audit$celltype)
  audit$n_control[a] <- n[1]
  audit$n_disease[a] <- n[2]
  if (any(n < 3L)) next
  x <- pb[, eligible$donor, drop = FALSE]
  model <- tryCatch(fit_ql(x, eligible$group), error = function(e) e)
  rows <- which(result$celltype == ct)
  if (inherits(model, "error")) {
    audit$status[a] <- result$status[rows] <- conditionMessage(model)
    next
  }
  models[[ct]] <- model
  audit$status[a] <- "Tested"
  t <- model$table[match(targets, rownames(model$table)), , drop = FALSE]
  result$log2FC[rows] <- t$logFC
  result$lnFC[rows] <- t$logFC * log(2)
  result$p_value[rows] <- t$PValue
  result$status[rows] <- ifelse(is.finite(t$PValue), "Tested", "Gene absent or filtered")
  for (omitted in eligible$donor) {
    keep <- eligible$donor != omitted
    n_sub <- table(factor(eligible$group[keep], levels = c("Normal", "DKD")))
    z <- data.frame(celltype = ct, gene = targets, omitted_donor = omitted,
      omitted_group = eligible$group[match(omitted, eligible$donor)],
      n_control = n_sub[1], n_disease = n_sub[2],
      lnFC = NA_real_, p_value = NA_real_, status = "Insufficient donors")
    if (all(n_sub >= 2L)) {
      refit <- tryCatch(fit_ql(x[, keep, drop = FALSE], eligible$group[keep],
                               retained = model$genes), error = function(e) e)
      if (inherits(refit, "error")) z$status <- conditionMessage(refit) else {
        t <- refit$table[match(targets, rownames(refit$table)), , drop = FALSE]
        z$lnFC <- t$logFC * log(2)
        z$p_value <- t$PValue
        z$status <- ifelse(is.finite(t$PValue), "Tested", "Gene absent or filtered")
      }
    }
    lodo[[length(lodo) + 1L]] <- z
  }
}
# Keep all 16 hypotheses even when a type is absent or cannot be modeled.
valid <- is.finite(result$p_value) & result$p_value >= 0 & result$p_value <= 1
result$FDR_16 <- p.adjust(ifelse(valid, result$p_value, 1), method = "BH")
result$FDR_16[!valid] <- NA_real_
stopifnot(nrow(result) == 16L)
write.csv(result, file.path(out, "target_results.csv"), row.names = FALSE)
write.csv(audit, file.path(out, "model_status.csv"), row.names = FALSE)
saveRDS(raw_sums, file.path(out, "raw_pseudobulk_counts.rds"))
saveRDS(models, file.path(out, "models_and_all_gene_tables.rds"))
lodo_table <- if (length(lodo)) do.call(rbind, lodo) else data.frame(
  celltype = character(), gene = character(), omitted_donor = character(),
  omitted_group = character(), n_control = integer(), n_disease = integer(),
  lnFC = numeric(), p_value = numeric(), status = character())
write.csv(lodo_table, file.path(out, "leave_one_donor_out.csv"), row.names = FALSE)

p <- ggplot(coverage, aes(donor, celltype, fill = log1p(n_cells))) +
  geom_tile(color = "white") + geom_text(aes(label = n_cells), size = 3) +
  scale_fill_gradient(low = "grey95", high = "#31688E") +
  labs(x = "Biological donor", y = NULL, fill = "log(1 + cells)") + theme_bw()
ggsave(file.path(out, "coverage.pdf"), p, width = 8, height = 5)
if (any(valid)) {
  p <- ggplot(result[valid, ], aes(celltype, lnFC, color = gene)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = position_dodge(width = 0.3), size = 3) +
    labs(x = NULL, y = "ln fold change: DKD versus Normal") + theme_bw()
  ggsave(file.path(out, "target_estimates.pdf"), p, width = 7, height = 4)
}
if (any(is.finite(lodo_table$lnFC))) {
  p <- ggplot(lodo_table, aes(omitted_donor, lnFC, color = omitted_group)) +
    geom_point(na.rm = TRUE) +
    geom_hline(data = result[valid, ], aes(yintercept = lnFC), linetype = 2) +
    facet_grid(gene ~ celltype, scales = "free_x") +
    labs(x = "Omitted donor", y = "Sensitivity ln fold change", color = "Group") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out, "leave_one_donor_out.pdf"), p, width = 10, height = 6)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
