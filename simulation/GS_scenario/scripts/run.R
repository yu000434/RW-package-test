#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run.R TASK_FILE TASK_ID RAW_DIR")

script <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
root <- normalizePath(file.path(dirname(script), "..", ".."))
core <- list.files(file.path(dirname(root), "R"), pattern = "\\.R$", full.names = TRUE)
invisible(lapply(core, source))
source(file.path(root, "results.R"))
source(file.path(root, "GS_scenario", "scripts", "generate_gs_data.R"))
source(file.path(root, "GS_scenario", "scripts", "run_one.R"))

task_id <- as.integer(args[[2L]])
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("TASK_ID must identify one task.")

raw <- do.call(rbind, lapply(seq_len(task$reps), function(i) {
  seed <- task$base_seed + i - 1L
  out <- run_one_gs(seed, task$n, task$validated, task$m, task$k, task$imputation)
  cbind(task[c("task_id", "run_id", "cell_id", "chunk_id")], out)
}))

dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
write.csv(raw, file.path(args[[3L]], sprintf("gs_task%04d.csv", task_id)),
          row.names = FALSE)
