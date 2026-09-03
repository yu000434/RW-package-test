#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript simulate_batch.R TASK_FILE TASK_ID RAW_DIR")
}
task_file <- args[[1L]]
task_id <- as.integer(args[[2L]])
raw_dir <- args[[3L]]
if (is.na(task_id)) stop("TASK_ID must be an integer.")

script_path <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())])
sim_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
source(file.path(sim_root, "GS_scenario", "PMM", "scripts", "run_one.R"))

tasks <- read.csv(task_file, stringsAsFactors = FALSE)
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L || task$scenario != "GS") stop("TASK_ID is not one GS task.")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(raw_dir, sprintf("gs_task%04d.csv", task_id))
if (file.exists(raw_path)) stop("Refusing to overwrite an existing task: ", raw_path)

seeds <- task$base_seed + seq_len(task$reps) - 1L
log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}
log_line("START run_id=", task$run_id, " task_id=", task_id,
         " n=", task$n, " validated=", task$validated, " m=", task$m,
         " k=", task$k, " reps=", task$reps,
         " seed_first=", min(seeds), " seed_last=", max(seeds))
raw <- do.call(rbind, lapply(seq_along(seeds), function(index) {
  seed <- seeds[[index]]
  log_line("PROGRESS task_id=", task_id, " replicate=", index, "/", task$reps,
           " seed=", seed, " state=START")
  out <- run_one_gs(
    seed = seed, n = task$n, validated = task$validated, m = task$m,
    k = task$k, sim_root = sim_root
  )
  log_line("PROGRESS task_id=", task_id, " replicate=", index, "/", task$reps,
           " seed=", seed, " state=DONE")
  cbind(task[rep(1L, nrow(out)), c("task_id", "run_id", "cell_id", "chunk_id")], out)
}))
temporary_path <- tempfile(pattern = "gs_task_", tmpdir = raw_dir)
write.csv(raw, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, raw_path)) {
  unlink(temporary_path)
  stop("Could not publish task output.")
}
log_line("COMPLETE run_id=", task$run_id, " task_id=", task_id,
         " raw_path=", normalizePath(raw_path))
