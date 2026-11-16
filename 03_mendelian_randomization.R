# Immune-cell cis-eQTL Mendelian randomization: a standardized example.
# Sources: OneK1K exposures; GCST90018832 outcome. Supply full input tables.
# The historical MR implementation was not recovered. This script illustrates
# the recorded thresholds; it does not certify the reported 126 associations.
# Input schemas and unresolved choices are documented in data/README.md.
suppressPackageStartupMessages(library(TwoSampleMR))
set.seed(123) # Example seed for stochastic sensitivity estimators.
out <- "results/mr"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Supply the actual ancestry-matched PLINK reference and executable.
ld_reference <- "data/mr/ld_reference"
plink_binary <- Sys.which("plink")
# Set the actual population prevalence used for liability-scale Steiger R2.
# This was not recorded; do not substitute the case fraction in the GWAS.
outcome_prevalence <- NA_real_
# Demonstration choice: described in the original text, but execution is unverified.
exclude_outcome_associated_snps <- TRUE

exposure <- read.csv("data/mr/exposure.csv", stringsAsFactors = FALSE)
outcome <- read.csv("data/mr/outcome.csv", stringsAsFactors = FALSE)
required <- c("SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval")
stopifnot(all(c(required, "gene", "cell_type", "samplesize") %in% names(exposure)),
          all(c(required, "ncase", "ncontrol") %in% names(outcome)),
          !anyDuplicated(outcome$SNP))
if (!is.finite(outcome_prevalence) || outcome_prevalence <= 0 ||
    outcome_prevalence >= 1) {
  stop("Set outcome_prevalence from the actual analysis specification.")
}
if (!nzchar(plink_binary) || !all(file.exists(paste0(ld_reference,
                                                       c(".bed", ".bim", ".fam"))))) {
  stop("Provide the PLINK executable and matching BED/BIM/FAM reference files.")
}
exposure$exposure_id <- paste(exposure$gene, exposure$cell_type, sep = "::")
exposure$phenotype <- exposure$exposure_id
lookup <- unique(exposure[c("exposure_id", "gene", "cell_type")])
stopifnot(!anyDuplicated(exposure[c("exposure_id", "SNP")] ))

# 1. Exposure association, LD and instrument-strength filters.
exposure <- subset(exposure, is.finite(beta) & is.finite(se) & se > 0 &
                     is.finite(pval) & pval < 1e-5 & samplesize > 2)
exp_dat <- format_data(exposure, type = "exposure", phenotype_col = "phenotype",
                       id_col = "exposure_id", samplesize_col = "samplesize")
exp_dat <- clump_data(exp_dat, clump_kb = 500, clump_r2 = 0.1, clump_p1 = 1,
                      bfile = ld_reference, plink_bin = plink_binary)
exp_dat$F_statistic <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
exp_dat <- subset(exp_dat, F_statistic > 10)
if (!nrow(exp_dat)) stop("No instruments passed exposure, LD and F filters.")
write.csv(exp_dat, file.path(out, "instruments_after_clumping.csv"), row.names = FALSE)

# 2. Harmonize alleles. action = 2 is an example palindromic-SNP setting.
outcome <- subset(outcome, SNP %in% exp_dat$SNP)
outcome$phenotype <- "DKD"
outcome$samplesize <- outcome$ncase + outcome$ncontrol
out_dat <- format_data(outcome, type = "outcome", phenotype_col = "phenotype",
                       samplesize_col = "samplesize", ncase_col = "ncase",
                       ncontrol_col = "ncontrol")
dat <- harmonise_data(exp_dat, out_dat, action = 2)
write.csv(dat, file.path(out, "harmonization_audit.csv"), row.names = FALSE)
dat <- subset(dat, mr_keep & is.finite(pval.outcome) &
                is.finite(eaf.outcome) & eaf.outcome > 0 & eaf.outcome < 1 &
                is.finite(se.outcome) & se.outcome > 0 &
                ncase.outcome > 0 & ncontrol.outcome > 0)
if (exclude_outcome_associated_snps) dat <- subset(dat, pval.outcome >= 5e-6)
if (!nrow(dat)) stop("No instruments remained after harmonization/outcome filtering.")

# 3. Steiger direction filtering; retain the full audit including test P values.
dat$rsq.exposure <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure,
                                   dat$samplesize.exposure)^2
dat$rsq.outcome <- get_r_from_lor(dat$beta.outcome, dat$eaf.outcome,
                                 dat$ncase.outcome, dat$ncontrol.outcome,
                                 prevalence = outcome_prevalence)^2
dat$effective_n.exposure <- dat$samplesize.exposure
dat$effective_n.outcome <- effective_n(dat$ncase.outcome, dat$ncontrol.outcome)
dat <- steiger_filtering(dat)
write.csv(dat, file.path(out, "steiger_all.csv"), row.names = FALSE)
dat <- subset(dat, !is.na(steiger_dir) & steiger_dir)
if (!nrow(dat)) stop("No instruments passed Steiger direction filtering.")
write.csv(dat, file.path(out, "harmonized_instruments.csv"), row.names = FALSE)

# 4. Wald ratio for one instrument; IVW for two or more instruments.
# Additional estimators are returned only where the method is estimable.
methods <- c("mr_wald_ratio", "mr_ivw", "mr_egger_regression",
             "mr_weighted_median", "mr_weighted_mode", "mr_simple_mode")
estimates <- mr(dat, method_list = methods)
estimates <- merge(estimates, lookup, by.x = "id.exposure", by.y = "exposure_id")
estimates$OR <- exp(estimates$b)
estimates$OR_low <- exp(estimates$b - 1.96 * estimates$se)
estimates$OR_high <- exp(estimates$b + 1.96 * estimates$se)
write.csv(estimates, file.path(out, "all_estimators.csv"), row.names = FALSE)
primary <- subset(estimates,
                  (nsnp == 1 & method == "Wald ratio") |
                    (nsnp >= 2 & method == "Inverse variance weighted"))
# Example multiplicity rule: BH across ALL estimable primary gene-cell tests.
# The original adjustment family is unverified. Never adjust only 126 selected rows.
primary$FDR <- p.adjust(primary$pval, method = "BH")
write.csv(primary, file.path(out, "primary_results.csv"), row.names = FALSE)
selected <- subset(primary, !is.na(FDR) & FDR < 0.05)
write.csv(selected, file.path(out, "selected_results_example.csv"), row.names = FALSE)
writeLines(sort(unique(selected$gene)), file.path(out, "prioritized_genes_example.txt"))

# 5. Sensitivity outputs. Missing estimates indicate an inapplicable method.
write.csv(mr_heterogeneity(dat, method_list = c("mr_ivw", "mr_egger_regression")),
          file.path(out, "heterogeneity.csv"), row.names = FALSE)
write.csv(mr_pleiotropy_test(dat), file.path(out, "egger_intercept.csv"), row.names = FALSE)
counts <- table(dat$id.exposure)
loo_dat <- subset(dat, id.exposure %in% names(counts[counts >= 3]))
if (nrow(loo_dat)) {
  write.csv(mr_leaveoneout(loo_dat, method = mr_ivw),
            file.path(out, "leave_one_variant_out.csv"), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
