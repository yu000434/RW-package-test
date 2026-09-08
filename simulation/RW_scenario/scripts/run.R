#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run.R TASK_FILE TASK_ID RAW_DIR")

script <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(script), "..", ".."))
core <- list.files(file.path(dirname(root), "R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(core, source))
source(file.path(root, "results.R"))
source(file.path(root, "RW_scenario", "scripts", "generate_rw_data.R"))
source(file.path(root, "RW_scenario", "scripts", "run_one.R"))

task_id <- as.integer(args[[2L]])
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("TASK_ID must identify one task.")

scenarios <- c(S1 = "robins_1", S2a = "robins_2_1", S2b = "robins_2_2",
               S2c = "robins_2_3", S3a = "robins_3_1", S3b = "robins_3_2")
raw <- do.call(rbind, lapply(seq_len(task$reps), function(i) {
  do.call(rbind, lapply(seq_along(scenarios), function(j) {
    offset <- if (task$imputation == "pmm") (j - 1L) * 100000L else 0L
    seed <- task$base_seed + i - 1L + offset
    out <- run_one_rw(seed, scenarios[[j]], task$n, task$m, task$k, task$imputation)
    cbind(task[c("task_id", "run_id", "chunk_id")],
          scenario_label = names(scenarios)[[j]], out)
  }))
}))

dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
write.csv(raw, file.path(args[[3L]], sprintf("rw_task%04d.csv", task_id)),
          row.names = FALSE)
