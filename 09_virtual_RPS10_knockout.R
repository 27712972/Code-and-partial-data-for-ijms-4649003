#!/usr/bin/env Rscript
# Illustrative implementation of the manuscript's human renal T-cell analysis.
# This is not a claim that the manuscript figures have been reproduced.
# Run from the repository root, or pass input_RDS and output_directory as arguments.
# Requires scTenifoldKnk >= 2.0.0 and scTenifoldNet >= 1.4.3.
# Network FDR and direction scores are model predictions, not measured expression.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[1L] else
  "results/single_cell_processing/human/human_kidney_processed.rds"
output_dir <- if (length(args) >= 2L) args[2L] else "results/virtual_knockout"
cfg <- list(target = "RPS10", groups = c("Normal", "DKD"),
  min_cells_per_donor = 20L, min_donors = 3L, max_cells_per_donor = 300L,
  min_detection = 0.05, n_hvg = 1500L, seeds = c(101L, 202L, 303L),
  n_networks = 10L, cells_per_network = 500L, n_components = 3L,
  tensor_rank = 3L, tensor_max_iter = 1000L, edge_quantile = 0.9,
  manifold_dim = 2L, n_cores = 1L, fdr_cutoff = 0.05)
packages <- c("Seurat", "SeuratObject", "Matrix", "scTenifoldNet", "scTenifoldKnk")
missing_packages <- packages[!vapply(packages, requireNamespace,
  logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ",
  paste(missing_packages, collapse = ", "))
for (p in c("scTenifoldKnk", "scTenifoldNet")) {
  required <- if (p == "scTenifoldKnk") "2.0.0" else "1.4.3"
  if (utils::compareVersion(as.character(utils::packageVersion(p)), required) < 0L)
    stop(p, " must be version ", required, " or newer.")
}
if (!"regressLibSize" %in% names(formals(scTenifoldKnk::knockoutDirection)))
  stop("The installed knockoutDirection API lacks regressLibSize.")
if (!file.exists(input_file)) stop("Input RDS not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
saveRDS(cfg, file.path(output_dir, "parameters.rds"))
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
is_ribosomal <- function(x) grepl("^(RP[LS][0-9]|RPLP[0-9]|RPSA$|MRP[LS][0-9])", x)
row_variance <- function(x) {
  n <- ncol(x)
  pmax(0, (Matrix::rowSums(x^2) - Matrix::rowSums(x)^2 / n) / (n - 1L))
}
sample_cells <- function(meta, seed) {
  set.seed(seed)
  unlist(lapply(split(meta$cell, meta$donor), function(x)
    sample(x, min(length(x), cfg$max_cells_per_donor), replace = FALSE)),
    use.names = FALSE)
}

# Extract all RNA genes for T cells; retain the full library for CPM denominators.
obj <- readRDS(input_file)
if (!inherits(obj, "Seurat") || !"RNA" %in% SeuratObject::Assays(obj))
  stop("Input must be a Seurat object with an RNA assay.")
meta <- obj[[]]
if (!all(c("celltype", "group", "donor") %in% names(meta)))
  stop("Metadata must contain celltype, group, and donor.")
keep <- !is.na(meta$celltype) & meta$celltype == "T cells" &
  !is.na(meta$group) & meta$group %in% cfg$groups
if (!any(keep)) stop("No T cells with the requested group labels were found.")
obj <- subset(obj, cells = rownames(meta)[keep])
meta <- obj[[]]
meta <- data.frame(cell = rownames(meta), donor = as.character(meta$donor),
  group = as.character(meta$group), stringsAsFactors = FALSE)
if (anyNA(meta) || any(!nzchar(meta$donor))) stop("Missing donor metadata.")
if (any(vapply(split(meta$group, meta$donor), function(x)
  length(unique(x)), integer(1)) != 1L)) stop("A donor belongs to multiple groups.")
rownames(meta) <- meta$cell
if (inherits(obj[["RNA"]], "Assay5")) {
  layers <- SeuratObject::Layers(obj[["RNA"]], search = "^counts($|\\.)")
  if (!length(layers)) stop("RNA raw counts are unavailable.")
  if (!identical(layers, "counts")) obj <- SeuratObject::JoinLayers(
    obj, assay = "RNA", layers = "counts", new = "counts")
  X <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
} else X <- SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts")
X <- methods::as(X[, meta$cell, drop = FALSE], "dgCMatrix")
if (anyDuplicated(rownames(X)) || anyDuplicated(colnames(X)))
  stop("Duplicate gene or cell identifiers.")
if (!length(X@x) || any(!is.finite(X@x)) || any(X@x < 0) ||
  any(abs(X@x - round(X@x)) > 1e-6) || any(Matrix::colSums(X) == 0))
  stop("Counts must be finite, nonnegative integers with nonzero cell libraries.")
if (!cfg$target %in% rownames(X)) stop("RPS10 is absent; use verified human gene symbols.")
coverage <- stats::aggregate(cell ~ donor + group, meta, length)
names(coverage)[3L] <- "n_cells"
coverage$eligible <- coverage$n_cells >= cfg$min_cells_per_donor
write_csv(coverage, file.path(output_dir, "donor_coverage.csv"))
meta <- meta[meta$donor %in% coverage$donor[coverage$eligible], , drop = FALSE]
for (g in cfg$groups) if (length(unique(meta$donor[meta$group == g])) < cfg$min_donors)
  stop("Fewer than three eligible donors in group ", g)
X <- X[, meta$cell, drop = FALSE]
write_csv(meta, file.path(output_dir, "included_Tcell_metadata.csv"))
sample_summary <- do.call(rbind, lapply(cfg$groups, function(g) {
  m <- meta[meta$group == g, , drop = FALSE]
  cells <- sample_cells(m, cfg$seeds[1L])
  if (length(cells) < cfg$cells_per_network)
    stop("Fewer than 500 donor-capped T cells in group ", g)
  available <- table(m$donor); retained <- table(m[cells, "donor"])
  data.frame(group = g, donor = names(available), available = as.integer(available),
    retained_seed101 = as.integer(retained[names(available)]))
}))
write_csv(sample_summary, file.path(output_dir, "network_sample_sizes.csv"))
cp <- scTenifoldNet::cpmNormalization(X)
expression <- data.frame(meta, RPS10_log1p_CPM = log1p(as.numeric(cp[cfg$target, ])))
write_csv(stats::aggregate(RPS10_log1p_CPM ~ donor + group, expression, mean),
  file.path(output_dir, "RPS10_donor_mean_logCPM.csv"))
write_csv(data.frame(group = cfg$groups, detection_fraction = vapply(cfg$groups,
  function(g) mean(X[cfg$target, meta$group == g] > 0), numeric(1))),
  file.path(output_dir, "RPS10_pooled_detection_before_subsampling.csv"))

# Select one common feature set, independently of knockout significance.
eligible <- rownames(X)
for (g in cfg$groups) {
  z <- X[, meta$cell[meta$group == g], drop = FALSE]
  good <- Matrix::rowMeans(z > 0) >= cfg$min_detection &
    row_variance(scTenifoldNet::cpmNormalization(z)) > .Machine$double.eps
  eligible <- intersect(eligible, rownames(z)[good])
}
eligible <- eligible[!grepl("^MT-", eligible)]
if (!cfg$target %in% eligible) stop("RPS10 did not pass detection/variance filters.")
feature_cells <- unlist(lapply(cfg$groups, function(g)
  sample_cells(meta[meta$group == g, , drop = FALSE], cfg$seeds[1L])), use.names = FALSE)
tmp <- SeuratObject::CreateSeuratObject(X[, feature_cells, drop = FALSE],
  min.cells = 0, min.features = 0)
tmp <- Seurat::FindVariableFeatures(tmp, selection.method = "vst",
  nfeatures = min(nrow(X), cfg$n_hvg * 3L), verbose = FALSE)
hvg <- head(intersect(SeuratObject::VariableFeatures(tmp), eligible), cfg$n_hvg)
genes <- unique(c(hvg, eligible[is_ribosomal(eligible)], cfg$target))
write_csv(data.frame(gene = genes, ribosomal = is_ribosomal(genes)),
  file.path(output_dir, "primary_network_features.csv"))

build_network <- function(counts, features, seed) {
  Z <- scTenifoldNet::cpmNormalization(counts)[features, , drop = FALSE]
  if (row_variance(Z[cfg$target, , drop = FALSE]) <= .Machine$double.eps)
    stop("RPS10 has insufficient variance in the sampled cells.")
  set.seed(seed)
  nets <- scTenifoldNet::makeNetworks(X = Z, nNet = cfg$n_networks,
    nCells = cfg$cells_per_network, nComp = cfg$n_components, scaleScores = TRUE,
    symmetric = FALSE, q = cfg$edge_quantile, nCores = cfg$n_cores)
  set.seed(seed)
  td <- scTenifoldNet::tensorDecomposition(xList = nets, K = cfg$tensor_rank,
    maxError = 1e-5, maxIter = cfg$tensor_max_iter, nDecimal = 3L, seed = seed)
  W <- as.matrix(td$X); diag(W) <- 0; W <- t(W)
  W <- W[features, features, drop = FALSE]
  if (any(!is.finite(W)) || sum(abs(W[cfg$target, ])) == 0)
    stop("The network is invalid or RPS10 has no outgoing edges.")
  W
}
run_ko <- function(W, counts, seed) {
  perturbed <- W; perturbed[cfg$target, ] <- 0
  set.seed(seed)
  M <- scTenifoldNet::manifoldAlignment(W, perturbed,
    d = cfg$manifold_dim, nCores = cfg$n_cores)
  D <- scTenifoldKnk::dRegulation(M, gKO = cfg$target, empiricalNull = FALSE)
  if (!all(c("gene", "distance", "p.adj") %in% names(D)) ||
    any(!is.finite(D$distance)) || any(!is.finite(D$p.adj)))
    stop("Invalid output from dRegulation.")
  names(D)[names(D) == "FC"] <- "network_relative_squared_distance"
  names(D)[names(D) == "p.value"] <- "network_p_value"
  names(D)[names(D) == "p.adj"] <- "network_FDR_BH"
  D$is_target <- D$gene == cfg$target
  D$above_numerical_noise <- D$distance > sqrt(.Machine$double.eps) * max(abs(M))
  for (adjust in c(FALSE, TRUE)) {
    scores <- scTenifoldKnk::knockoutDirection(X = counts, gKO = cfg$target,
      genes = rownames(W), t = 5, regressLibSize = adjust)
    s <- as.numeric(scores[1L, match(D$gene, colnames(scores))]); s[D$is_target] <- NA_real_
    suffix <- if (adjust) "depth_adjusted" else "unadjusted"
    D[[paste0("direction_score_", suffix)]] <- s
    D[[paste0("predicted_direction_", suffix)]] <- ifelse(!is.finite(s) | s == 0,
      "undetermined", ifelse(s < 0, "down", "up"))
  }
  D$selected <- !D$is_target & D$above_numerical_noise &
    is.finite(D$network_FDR_BH) & D$network_FDR_BH < cfg$fdr_cutoff
  D[order(D$network_FDR_BH, -D$distance), , drop = FALSE]
}

# Repeat primary networks; rebuild the ribosomal-exclusion network only at seed 101.
for (g in cfg$groups) {
  m <- meta[meta$group == g, , drop = FALSE]
  group_dir <- file.path(output_dir, g); dir.create(group_dir, showWarnings = FALSE)
  for (seed in cfg$seeds) {
    run_dir <- file.path(group_dir, paste0("seed_", seed))
    dir.create(run_dir, showWarnings = FALSE)
    cells <- sample_cells(m, seed); counts <- X[, cells, drop = FALSE]
    write_csv(m[cells, , drop = FALSE], file.path(run_dir, "selected_cells.csv"))
    W <- build_network(counts, genes, seed); D <- run_ko(W, counts, seed)
    saveRDS(W, file.path(run_dir, "primary_network.rds"))
    write_csv(D, file.path(run_dir, "primary_KO_all_genes.csv"))
    if (seed != cfg$seeds[1L]) next
    restricted <- genes[!is_ribosomal(genes) | genes == cfg$target]
    W2 <- build_network(counts, restricted, seed); D2 <- run_ko(W2, counts, seed)
    saveRDS(W2, file.path(run_dir, "ribosomal_exclusion_network.rds"))
    write_csv(D2, file.path(group_dir, "remove_other_ribosomal_KO.csv"))
    shared <- setdiff(intersect(D$gene, D2$gene), cfg$target)
    if (length(shared) < 3L) stop("Too few shared downstream genes for rank comparison.")
    comparison <- data.frame(gene = shared,
      primary_distance = D$distance[match(shared, D$gene)],
      exclusion_distance = D2$distance[match(shared, D2$gene)],
      predicted_direction = D2$predicted_direction_depth_adjusted[match(shared, D2$gene)])
    comparison$primary_rank_fraction <- rank(comparison$primary_distance) / nrow(comparison)
    comparison$exclusion_rank_fraction <- rank(comparison$exclusion_distance) / nrow(comparison)
    k <- min(100L, nrow(comparison))
    top_primary <- shared[head(order(-comparison$primary_distance, shared), k)]
    top_exclusion <- shared[head(order(-comparison$exclusion_distance, shared), k)]
    comparison$top100_primary <- shared %in% top_primary
    comparison$top100_exclusion <- shared %in% top_exclusion
    write_csv(comparison, file.path(group_dir, "shared_gene_rank_comparison.csv"))
    write_csv(data.frame(n_shared = length(shared), top_k = k,
      spearman = stats::cor(comparison$primary_distance, comparison$exclusion_distance,
        method = "spearman"), overlap = length(intersect(top_primary, top_exclusion)),
      jaccard = length(intersect(top_primary, top_exclusion)) /
        length(union(top_primary, top_exclusion))),
      file.path(group_dir, "ribosomal_sensitivity_rank_agreement.csv"))
    if (g == "DKD") {
      down <- D2[D2$selected & D2$predicted_direction_depth_adjusted == "down", , drop = FALSE]
      write_csv(down, file.path(output_dir, "ko_down_gene_details.csv"))
      writeLines(down$gene, file.path(output_dir, "ko_down_genes.txt"))
      writeLines(shared, file.path(output_dir, "tested_non_target_genes.txt"))
    }
  }
}
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Virtual knockout analysis finished. Inspect actual outputs and model assumptions.")
