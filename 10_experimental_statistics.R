# Illustrative analysis of animal-level validation measurements (Figure 4).
# Run from the package root: Rscript scripts/10_experimental_statistics.R
# Required package: ggplot2. No example observations are generated.
# Input: data/experiments/animal_measurements.csv
# Columns: animal_id, group, endpoint, value; group must be Con or Mod.
# Each endpoint must contain one value per animal from five animals per group.
# Include measurement units in endpoint names, for example UACR_mg_per_g.
# Values must already reflect the documented animal-level normalization.
# Average technical replicates within animals before creating this input.
# Raw IF/WB images and their quantification cannot be reconstructed here.
# For qPCR, confirm whether the original test used delta-Ct or relative
# expression; this script tests the supplied scale and does not choose it.
# These t-tests illustrate the stated method, not the original Prism settings.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 before running this script.")
}
input_file <- "data/experiments/animal_measurements.csv"
output_dir <- "results/experiments"
if (!file.exists(input_file)) stop("Missing input: ", input_file)
dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("animal_id", "group", "endpoint", "value")
if (!all(required %in% names(dat))) {
  stop("Input must contain: ", paste(required, collapse = ", "))
}
dat <- dat[, required, drop = FALSE]
if (nrow(dat) == 0L || anyNA(dat)) stop("Input is empty or contains missing values.")
dat$animal_id <- trimws(as.character(dat$animal_id))
dat$group <- trimws(as.character(dat$group))
dat$endpoint <- trimws(as.character(dat$endpoint))
dat$value <- suppressWarnings(as.numeric(dat$value))
if (any(!is.finite(dat$value))) stop("Every value must be numeric and finite.")
if (any(!nzchar(dat$animal_id)) || any(!nzchar(dat$endpoint))) {
  stop("Animal IDs and endpoint names must not be blank.")
}
if (!setequal(unique(dat$group), c("Con", "Mod"))) {
  stop("Groups must be exactly Con and Mod.")
}
if (anyDuplicated(dat[, c("animal_id", "endpoint")])) {
  stop("Keep one normalized value per animal per endpoint.")
}
assignment <- split(dat$group, dat$animal_id)
if (any(vapply(assignment, function(x) length(unique(x)) != 1L, logical(1)))) {
  stop("An animal ID must belong to only one group.")
}
animals <- unique(dat[, c("animal_id", "group")])
if (any(table(factor(animals$group, levels = c("Con", "Mod"))) != 5L)) {
  stop("This example expects five distinct validation animals per group.")
}
dat$group <- factor(dat$group, levels = c("Con", "Mod"))
counts <- table(dat$endpoint, dat$group)
if (any(counts != 5L)) stop("Each endpoint requires five animal values per group.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summary_rows <- list()
test_rows <- list()
endpoint_names <- sort(unique(dat$endpoint))
for (i in seq_along(endpoint_names)) {
  endpoint_name <- endpoint_names[i]
  d <- dat[dat$endpoint == endpoint_name, , drop = FALSE]
  s <- do.call(rbind, lapply(levels(d$group), function(g) {
    v <- d$value[d$group == g]
    data.frame(endpoint = endpoint_name, group = g, n = length(v),
               mean = mean(v), sd = sd(v))
  }))
  s$group <- factor(s$group, levels = c("Con", "Mod"))
  mod <- d$value[d$group == "Mod"]
  con <- d$value[d$group == "Con"]
  fit <- tryCatch(t.test(mod, con, paired = FALSE, var.equal = TRUE,
                        alternative = "two.sided", conf.level = 0.95),
                  error = function(e) stop(endpoint_name, ": ", conditionMessage(e)))
  summary_rows[[i]] <- s
  test_rows[[i]] <- data.frame(
    endpoint = endpoint_name, contrast = "Mod - Con", n_Con = length(con),
    n_Mod = length(mod), mean_difference = mean(mod) - mean(con),
    ci_low = unname(fit$conf.int[1]), ci_high = unname(fit$conf.int[2]),
    t_statistic = unname(fit$statistic), df = unname(fit$parameter),
    p_value = fit$p.value, method = "Two-sided unpaired Student t-test",
    multiple_testing_adjustment = "None")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = value, color = group)) +
    ggplot2::geom_point(size = 2.5,
      position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 101)) +
    ggplot2::geom_errorbar(data = s,
      ggplot2::aes(x = group, ymin = mean - sd, ymax = mean + sd),
      inherit.aes = FALSE, width = 0.20, color = "black", linewidth = 0.6) +
    ggplot2::geom_point(data = s, ggplot2::aes(x = group, y = mean),
      inherit.aes = FALSE, shape = 18, size = 3, color = "black") +
    ggplot2::scale_color_manual(values = c(Con = "#3977B7", Mod = "#C64E52")) +
    ggplot2::labs(x = NULL, y = endpoint_name, title = endpoint_name,
      subtitle = paste0("Unpaired t-test, P = ", format.pval(fit$p.value, digits = 3)),
      caption = "Individual animals; black diamond and error bars: mean +/- SD") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", endpoint_name)
  figure_file <- sprintf("%02d_%s.pdf", i, safe_name)
  ggplot2::ggsave(file.path(output_dir, figure_file), p, width = 4.5, height = 4.3)
}
write.csv(do.call(rbind, summary_rows), file.path(output_dir, "summary.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, test_rows), file.path(output_dir, "t_tests.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Saved animal-level summaries, t-tests and endpoint figures to ", output_dir)
