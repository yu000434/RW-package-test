#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript summarize.R TASK_FILE RAW_DIR SUMMARY_FILE")
}
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
raw_dir <- args[[2L]]
summary_file <- args[[3L]]
paths <- file.path(raw_dir, sprintf("rw_task%04d.csv", tasks$task_id))
if (any(!file.exists(paths))) stop("All task files are required before summarizing.")
raw <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
if (nrow(raw) != 6L * sum(tasks$reps) ||
    anyDuplicated(raw[c("scenario", "seed")])) {
  stop("The raw RW results are incomplete or duplicated.")
}
if (any(!is.finite(raw$estimate)) || any(raw$rb_total_var <= 0) ||
    any(raw$rr_total_var <= 0)) {
  stop("The raw RW results contain an invalid estimate or variance.")
}

summarize_cell <- function(x) {
  if (nrow(x) < 2L) stop("At least two replications are required for a summary.")
  center <- mean(x$estimate)
  coverage_summary <- function(hit) {
    interval <- stats::binom.test(sum(hit), length(hit))$conf.int
    c(estimate = mean(hit), lower = interval[[1L]], upper = interval[[2L]])
  }
  rb_empirical <- coverage_summary(abs(x$estimate - center) <= 1.96 * x$rb_se)
  rr_empirical <- coverage_summary(abs(x$estimate - center) <= 1.96 * x$rr_se)
  rb_beta0 <- coverage_summary(abs(x$estimate - x$beta0) <= 1.96 * x$rb_se)
  rr_beta0 <- coverage_summary(abs(x$estimate - x$beta0) <= 1.96 * x$rr_se)
  data.frame(
    scenario = x$scenario[[1L]], n = x$n[[1L]], validated = NA_integer_,
    m = x$m[[1L]], k = x$k[[1L]], reps = nrow(x),
    mean_estimate = center, bias = mean(x$estimate - x$beta0), empirical_sd = stats::sd(x$estimate),
    mean_rb_se = mean(x$rb_se), mean_rr_se = mean(x$rr_se),
    rb_se_empirical_sd = mean(x$rb_se) / stats::sd(x$estimate),
    rr_se_empirical_sd = mean(x$rr_se) / stats::sd(x$estimate),
    rb_coverage_empirical_center = rb_empirical[["estimate"]],
    rb_coverage_empirical_center_lower = rb_empirical[["lower"]],
    rb_coverage_empirical_center_upper = rb_empirical[["upper"]],
    rr_coverage_empirical_center = rr_empirical[["estimate"]],
    rr_coverage_empirical_center_lower = rr_empirical[["lower"]],
    rr_coverage_empirical_center_upper = rr_empirical[["upper"]],
    rb_coverage_beta0 = rb_beta0[["estimate"]],
    rb_coverage_beta0_lower = rb_beta0[["lower"]],
    rb_coverage_beta0_upper = rb_beta0[["upper"]],
    rr_coverage_beta0 = rr_beta0[["estimate"]],
    rr_coverage_beta0_lower = rr_beta0[["lower"]],
    rr_coverage_beta0_upper = rr_beta0[["upper"]],
    mean_omega_var = mean(x$omega_var), mean_kappa_var = mean(x$kappa_var),
    mean_cross_var = mean(x$cross_var), mean_kappa_norm = mean(x$kappa_norm),
    mean_u_bar_norm = mean(x$u_bar_norm), mean_s_mis_imp_norm = mean(x$s_mis_imp_norm),
    mean_max_donor_reuse = mean(x$max_donor_reuse)
  )
}

summary <- do.call(rbind, lapply(split(raw, raw$scenario), summarize_cell))
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, summary_file, row.names = FALSE)
cat("SUMMARY_PATH=", normalizePath(summary_file), "\n", sep = "")
