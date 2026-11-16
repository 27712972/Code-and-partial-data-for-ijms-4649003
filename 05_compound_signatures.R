# Illustrative Enrichr query; run from the repository root. Internet access required.
# Signature up/down refers to the reference transcriptional response to exposure.
# Enrichment does not establish direct compound-target binding or treatment efficacy.
# Official API: https://maayanlab.cloud/Enrichr/templates/help/api-content.html
required <- c("httr", "jsonlite", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors = FALSE)
input <- "data/mr/prioritized_genes.txt"
selection_file <- "data/compound_signatures/selected_terms.csv"
out <- "results/compound_signatures"
base_url <- "https://maayanlab.cloud/Enrichr"
libraries <- c("NIBR_DRUGseq_2025_up", "NIBR_DRUGseq_2025_down")
if (!file.exists(input)) stop("Supply the exact retained MR gene list: ", input)
genes <- unique(trimws(readLines(input, warn = FALSE)))
genes <- genes[nzchar(genes)]
if (!length(genes)) stop("The MR input gene list is empty.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(genes, file.path(out, "submitted_symbols.txt"))

save_response <- function(response, name) {
  raw <- httr::content(response, as = "text", encoding = "UTF-8")
  writeLines(raw, file.path(out, name), useBytes = TRUE)
  httr::stop_for_status(response)
  jsonlite::fromJSON(raw, simplifyVector = FALSE)
}
metadata <- save_response(httr::GET(paste0(base_url, "/datasetStatistics"), httr::timeout(120)),
                          "library_metadata.json")
available <- vapply(metadata$statistics, function(x) x$libraryName, character(1))
if (!all(libraries %in% available)) stop("Required NIBR libraries are unavailable in Enrichr.")
# Running this POST submits the supplied gene symbols to the public Enrichr service.
submission <- save_response(httr::POST(paste0(base_url, "/addList"), encode = "multipart",
  body = list(list = paste(genes, collapse = "\n"), description = "IJMS MR gene set"),
  httr::timeout(120)), "submission.json")
if (is.null(submission$userListId)) stop("Enrichr did not return a userListId.")
empty <- data.frame(library = character(), term = character(), rank = integer(),
  p_value = numeric(), adjusted_p_value = numeric(), combined_score = numeric(),
  overlap_count = integer(), overlapping_genes = character(), direction = character())
tables <- list()
for (lib in libraries) {
  response <- save_response(httr::GET(paste0(base_url, "/enrich"),
    query = list(userListId = submission$userListId, backgroundType = lib), httr::timeout(120)),
    paste0(lib, "_full_response.json"))
  if (!(lib %in% names(response))) stop("Unexpected response schema for ", lib)
  rows <- response[[lib]]
  table <- if (!length(rows)) empty else do.call(rbind, lapply(rows, function(row) {
    if (length(row) < 7L) stop("Unexpected enrichment row schema for ", lib)
    overlaps <- as.character(unlist(row[[6]], use.names = FALSE))
    data.frame(library = lib, term = row[[2]], rank = as.integer(row[[1]]),
      p_value = as.numeric(row[[3]]), adjusted_p_value = as.numeric(row[[7]]),
      combined_score = as.numeric(row[[5]]), overlap_count = length(overlaps),
      overlapping_genes = paste(overlaps, collapse = ";"),
      direction = if (grepl("_up$", lib)) "Reference upregulated" else "Reference downregulated")
  }))
  tables[[lib]] <- table
  write.csv(table, file.path(out, paste0(lib, "_all_returned.csv")), row.names = FALSE)
}
all_results <- do.call(rbind, tables)
# Keep Enrichr's adjusted P values; do not recompute BH after filtering the response.
significant <- all_results[!is.na(all_results$adjusted_p_value) &
                           all_results$adjusted_p_value < 0.05, ]
write.csv(significant, file.path(out, "all_significant_signatures.csv"), row.names = FALSE)
writeLines(c(paste("Retrieved at UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("Libraries:", paste(libraries, collapse = ", ")),
  "Live service results may differ from archived original analyses.",
  "No automatic selection of the seven displayed manuscript signatures is implemented.",
  "Selected signatures require exact library and term values in selected_terms.csv.",
  "Up/down denotes a reference transcriptional response, not direct molecular binding.",
  capture.output(sessionInfo())), file.path(out, "run_notes.txt"))

if (file.exists(selection_file)) {
  selected <- read.csv(selection_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("library", "term") %in% names(selected)) || !nrow(selected)) {
    stop("selected_terms.csv must contain nonempty library and term columns.")
  }
  key <- function(d) paste(d$library, d$term, sep = "\r")
  if (anyDuplicated(key(selected))) stop("Duplicate selected library/term pairs.")
  index <- match(key(selected), key(all_results))
  if (anyNA(index)) stop("Selected signatures absent from response: ",
                         paste(selected$term[is.na(index)], collapse = "; "))
  display <- all_results[index, , drop = FALSE]
  write.csv(display, file.path(out, "manually_selected_signatures.csv"), row.names = FALSE)
  display$label <- paste(display$term, display$direction, sep = " | ")
  display$label <- factor(display$label,
                          levels = display$label[order(display$combined_score)])
  p <- ggplot(display, aes(combined_score, label, fill = direction)) + geom_col() +
    labs(x = "Enrichr combined score", y = NULL, fill = "Reference response") +
    theme_bw(base_size = 10)
  ggsave(file.path(out, "manually_selected_signatures.pdf"), p,
         width = 12, height = max(3, nrow(display) * 0.45 + 1.5), limitsize = FALSE)
} else {
  message("No selected_terms.csv supplied; all significant results were exported without selection.")
}
