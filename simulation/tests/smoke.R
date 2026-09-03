#!/usr/bin/env Rscript

args <- commandArgs()
script <- sub("^--file=", "", args[grep("^--file=", args)])
sim <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)

compare_row <- function(out, path, seed, scenario = NULL) {
  ref <- read.csv(path, stringsAsFactors = FALSE)
  ref <- ref[ref$seed == seed, , drop = FALSE]
  if (!is.null(scenario)) ref <- ref[ref$scenario == scenario, , drop = FALSE]
  ref <- ref[, names(out), drop = FALSE]
  stopifnot(nrow(ref) == 1L)
  for (name in names(out)) {
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

method_files <- c("pmmrw.R", "topk.R", "rw_variance.R", "metrics.R", "rb.R")
invisible(lapply(file.path(sim, "R", method_files), source))
source(file.path(sim, "GS_scenario", "PMM", "scripts", "generate_gs_data.R"))
source(file.path(sim, "GS_scenario", "PMM", "scripts", "run_one.R"))
gs <- run_one_gs(8500000L, 2000L, 500L, 5L, 5L)
compare_row(gs, file.path(sim, "results", "gs5_5", "raw", "gs_task0001.csv"), 8500000L)
cat("GS_SMOKE=PASS\n")

source(file.path(sim, "RW_scenario", "PMM", "scripts", "generate_rw_data.R"))
source(file.path(sim, "RW_scenario", "PMM", "scripts", "run_one.R"))
rw <- run_one_rw(12000000L, "robins_1", 150L, 20L, 5L)
compare_row(rw, file.path(sim, "results", "rw150_20_5", "raw", "rw_task0001.csv"),
            12000000L, "robins_1")
cat("RW_SMOKE=PASS\n")
