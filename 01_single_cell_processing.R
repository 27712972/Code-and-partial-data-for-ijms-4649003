# Single-cell processing: mouse and human kidney
# Compact demonstration for the IJMS revision; not a verified reproduction.
# Reference analysis: Seurat 5.0.2. Required packages: Seurat, harmony, ggplot2.
# Run from the project root after editing the sample IDs and input paths.
# Input: one 10x gene-expression matrix directory per biological sample,
# with gene symbols in the second column of features.tsv(.gz).
# Mouse alignment (Cell Ranger 6.0, GRCm39) is upstream of this script.
# Doublet detection is not implemented because its original settings are absent.
# The QC cutoffs below follow the manuscript; human resolution 1.0 follows
# the revised analysis (RNA_snn_res.1), replacing the original value of 0.6.
# Example choices not fully specified in the manuscript: scale factor 10000,
# 2000 variable genes, seed 123, per-sample 97th percentile before QC,
# and Harmony with 20 dimensions for the human dataset.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

set.seed(123)
output_root <- "results/single_cell_processing"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Mouse sample IDs are placeholders; replace them with the actual IDs.
mouse_samples <- data.frame(
  sample = c("Con01", "Con02", "Con03", "Mod01", "Mod02", "Mod03"),
  group = rep(c("Con", "Mod"), each = 3)
)
human_samples <- data.frame(
  sample = c("NM01", "NM02", "NM03", "DKD01", "DKD02", "DKD03"),
  group = rep(c("Normal", "DKD"), each = 3)
)
mouse_samples$path <- file.path("data", "mouse", mouse_samples$sample)
human_samples$path <- file.path("data", "human", human_samples$sample)

process_kidney <- function(samples, species, resolution) {
  out <- file.path(output_root, species)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  qc_summary <- list()

  # 1. Read samples and apply species-specific quality control.
  objects <- lapply(seq_len(nrow(samples)), function(i) {
    counts <- Read10X(data.dir = samples$path[i], gene.column = 2)
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
    if (is.null(counts)) stop("A gene-expression count matrix is required.")
    obj <- CreateSeuratObject(
      counts = counts, project = samples$sample[i],
      min.cells = 0, min.features = 0
    )
    obj$sample_id <- samples$sample[i]
    obj$donor <- samples$sample[i]
    obj$group <- samples$group[i]
    mt_pattern <- if (species == "mouse") "^mt-" else "^MT-"
    if (!any(grepl(mt_pattern, rownames(obj)))) {
      stop("Mitochondrial genes were not found; check gene-symbol annotation.")
    }
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
    count_cutoff <- NA_real_

    if (species == "mouse") {
      keep <- obj$nFeature_RNA >= 200 & obj$percent.mt <= 25
    } else {
      hb_genes <- intersect(
        c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2",
          "HBM", "HBQ1", "HBZ"), rownames(obj)
      )
      if (!length(hb_genes)) {
        stop("Hemoglobin genes were not found; check gene-symbol annotation.")
      }
      obj[["percent.hb"]] <- PercentageFeatureSet(obj, features = hb_genes)
      count_cutoff <- unname(quantile(obj$nCount_RNA, probs = 0.97))
      keep <- obj$nFeature_RNA > 300 & obj$nFeature_RNA < 7000 &
        obj$nCount_RNA > 1000 & obj$nCount_RNA < count_cutoff &
        obj$percent.mt < 10 & obj$percent.hb < 3
    }
    keep[is.na(keep)] <- FALSE
    qc_summary[[i]] <<- data.frame(
      sample = samples$sample[i], group = samples$group[i],
      cells_before = ncol(obj), cells_after = sum(keep),
      count_97th_percentile = count_cutoff
    )
    if (!any(keep)) stop(paste("No cells passed QC:", samples$sample[i]))
    subset(obj, cells = colnames(obj)[keep])
  })
  write.csv(do.call(rbind, qc_summary), file.path(out, "qc_summary.csv"),
            row.names = FALSE)

  # 2. Combine samples within each species, retaining raw RNA counts.
  obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample)
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj, assay = "RNA")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 20,
                seed.use = 123, verbose = FALSE)

  # 3. Correct sample effects in the embedding and cluster cells.
  # Do not use the disease group as the batch variable.
  set.seed(123)
  obj <- RunHarmony(obj, group.by.vars = "sample_id", dims.use = 1:20,
                    verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20,
                       verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, random.seed = 123,
                      verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20,
                 seed.use = 123, verbose = FALSE)
  obj <- RunTSNE(obj, reduction = "harmony", dims = 1:20,
                 seed.use = 123)

  # 4. Apply curated labels only after checking markers for this exact run.
  # Optional input: data/<species>/cluster_annotation.csv; columns:
  # cluster,celltype. Use the harmonized markers in Supplementary Table S1.
  # If absent, export a template and retain cluster labels for plotting.
  clusters <- sort(unique(as.integer(as.character(obj$seurat_clusters))))
  annotation_file <- file.path("data", species, "cluster_annotation.csv")
  plot_group <- "seurat_clusters"
  obj$celltype <- "Unassigned"
  if (file.exists(annotation_file)) {
    annotation <- read.csv(annotation_file, colClasses = "character")
    stopifnot(all(c("cluster", "celltype") %in% names(annotation)),
              !anyDuplicated(annotation$cluster),
              !anyNA(annotation$celltype),
              all(nzchar(trimws(annotation$celltype))))
    index <- match(as.character(obj$seurat_clusters), annotation$cluster)
    if (anyNA(index)) stop("The annotation table must cover every cluster.")
    obj$celltype <- annotation$celltype[index]
    plot_group <- "celltype"
  } else {
    write.csv(data.frame(cluster = clusters, celltype = "Unassigned"),
              file.path(out, "cluster_annotation_template.csv"),
              row.names = FALSE)
  }

  # 5. Save the processed object, per-cell metadata and overview plots.
  # Cluster numbers are outputs; do not force the manuscript's cluster count.
  for (reduction in c("umap", "tsne")) {
    p <- DimPlot(obj, reduction = reduction, group.by = plot_group,
                 label = TRUE, repel = TRUE, pt.size = 0.1)
    ggsave(file.path(out, paste0(reduction, ".pdf")), p, width = 8, height = 6)
  }
  metadata <- data.frame(cell_id = colnames(obj), obj[[]], row.names = NULL)
  write.csv(metadata, file.path(out, "cell_metadata.csv"), row.names = FALSE)
  saveRDS(obj, file.path(out, paste0(species, "_kidney_processed.rds")))
  obj
}

# Process species independently; do not integrate mouse and human cells.
mouse <- process_kidney(mouse_samples, species = "mouse", resolution = 0.8)
human <- process_kidney(human_samples, species = "human", resolution = 1.0)
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
