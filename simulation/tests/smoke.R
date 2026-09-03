#!/usr/bin/env Rscript

args <- commandArgs()
script <- sub("^--file=", "", args[grep("^--file=", args)])
sim <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

compare_row <- function(out, path, seed, scenario = NULL) {
  ref <- read.csv(path, stringsAsFactors = FALSE)
  ref <- ref[ref$seed == seed, , drop = FALSE]
  if (!is.null(scenario)) ref <- ref[ref$scenario == scenario, , drop = FALSE]
  stopifnot(nrow(ref) == 1L)
  fields <- intersect(names(out), names(ref))
  for (name in fields) {
    if (is.numeric(out[[name]])) {
      stopifnot(identical(is.na(out[[name]]), is.na(ref[[name]])))
      keep <- !is.na(out[[name]])
      if (any(keep)) {
        stopifnot(isTRUE(all.equal(out[[name]][keep], ref[[name]][keep], tolerance = 1e-10)))
      }
    } else {
      stopifnot(identical(out[[name]], ref[[name]]))
    }
  }
}

method_files <- c("pmmrw.R", "pmm_score.R", "variance.R", "results.R")
invisible(lapply(file.path(sim, "R", method_files), source))
source(file.path(sim, "GS_scenario", "scripts", "generate_gs_data.R"))
source(file.path(sim, "GS_scenario", "scripts", "run_one.R"))
gs <- run_one_gs(8500000L, 2000L, 500L, 5L, 5L)
compare_row(gs, file.path(sim, "results", "gs5_5", "raw", "gs_task0001.csv"), 8500000L)
gs_parametric <- run_one_gs(303100000L, 2000L, 500L, 5L, imputation = "parametric")
gs_expected <- c(estimate = 0.768114092516404, rb_se = 0.111246300044035,
                 rr_se = 0.194773158820324, beta0 = 0.783182537588776,
                 mean_analysis_n = 537.8,
                 median_abs_score_component_corr = 0.968917291713887)
stopifnot(isTRUE(all.equal(unlist(gs_parametric[names(gs_expected)]), gs_expected,
                           tolerance = 1e-12, check.attributes = FALSE)))
cat("GS_PMM_AND_PARAMETRIC_SMOKE=PASS\n")

source(file.path(sim, "RW_scenario", "scripts", "generate_rw_data.R"))
source(file.path(sim, "RW_scenario", "scripts", "run_one.R"))
rw <- run_one_rw(12000000L, "robins_1", 150L, 20L, 5L)
compare_row(rw, file.path(sim, "results", "rw150_20_5", "raw", "rw_task0001.csv"),
            12000000L, "robins_1")
rw_parametric <- run_one_rw(1001L, "robins_1", 150L, 20L, imputation = "parametric")
rw_expected <- c(estimate = 0.79821491158947, rb_se = 0.131276675826991,
                 rr_se = 0.139343439579577)
stopifnot(isTRUE(all.equal(unlist(rw_parametric[names(rw_expected)]), rw_expected,
                           tolerance = 1e-12, check.attributes = FALSE)))
cat("RW_PMM_AND_PARAMETRIC_SMOKE=PASS\n")
