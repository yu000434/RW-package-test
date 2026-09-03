#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Usage: Rscript summarize.R TASK_FILE RAW_DIR SUMMARY_FILE [CELL_ID]")
}
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
raw_dir <- args[[2L]]
summary_file <- args[[3L]]
if (length(args) == 4L) {
  cell_id <- as.integer(args[[4L]])
  tasks <- tasks[tasks$cell_id == cell_id, , drop = FALSE]
}
paths <- file.path(raw_dir, sprintf("gs_task%04d.csv", tasks$task_id))
if (any(!file.exists(paths))) stop("All task files are required before summarizing.")
raw <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
if (nrow(raw) != sum(tasks$reps)) stop("The raw results are incomplete.")

summarize_cell <- function(x) {
  center <- mean(x$estimate)
  empirical_sd <- stats::sd(x$estimate)
  coverage_summary <- function(hit) {
    interval <- stats::binom.test(sum(hit), length(hit))$conf.int
    c(estimate = mean(hit), lower = interval[[1L]], upper = interval[[2L]])
  }
  rb_empirical <- coverage_summary(abs(x$estimate - center) <= 1.96 * x$rb_se)
  rr_empirical <- coverage_summary(abs(x$estimate - center) <= 1.96 * x$rr_se)
  rb_beta0 <- coverage_summary(abs(x$estimate - x$beta0) <= 1.96 * x$rb_se)
  rr_beta0 <- coverage_summary(abs(x$estimate - x$beta0) <= 1.96 * x$rr_se)
  data.frame(
    scenario = "GS", n = x$n[[1L]], validated = x$validated[[1L]],
    m = x$m[[1L]], k = x$k[[1L]], reps = nrow(x),
    mean_estimate = center, bias = mean(x$estimate - x$beta0), empirical_sd = empirical_sd,
    mean_rb_se = mean(x$rb_se), mean_rr_se = mean(x$rr_se),
    rb_se_empirical_sd = mean(x$rb_se) / empirical_sd,
    rr_se_empirical_sd = mean(x$rr_se) / empirical_sd,
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

key <- interaction(raw$n, raw$validated, raw$m, raw$k, drop = TRUE)
summary <- do.call(rbind, lapply(split(raw, key), summarize_cell))
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, summary_file, row.names = FALSE)
