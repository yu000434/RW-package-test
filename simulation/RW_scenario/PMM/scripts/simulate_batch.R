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
source(file.path(sim_root, "RW_scenario", "PMM", "scripts", "run_one.R"))
scenarios <- c(S1 = "robins_1", S2a = "robins_2_1", S2b = "robins_2_2",
               S2c = "robins_2_3", S3a = "robins_3_1", S3b = "robins_3_2")
seed_block <- 100000L

tasks <- read.csv(task_file, stringsAsFactors = FALSE)
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("TASK_ID is not one RW task.")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(raw_dir, sprintf("rw_task%04d.csv", task_id))
if (file.exists(raw_path)) stop("Refusing to overwrite an existing task: ", raw_path)

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  flush.console()
}
log_line("START run_id=", task$run_id, " task_id=", task_id,
         " n=", task$n, " m=", task$m,
         " k=", task$k, " reps=", task$reps,
         " seed_first=", task$base_seed,
         " seed_last=", task$base_seed + task$reps - 1L)
raw <- do.call(rbind, lapply(seq_len(task$reps), function(replicate_index) {
  do.call(rbind, lapply(seq_along(scenarios), function(scenario_index) {
    scenario_label <- names(scenarios)[[scenario_index]]
    scenario <- scenarios[[scenario_index]]
    seed <- task$base_seed + replicate_index - 1L + (scenario_index - 1L) * seed_block
    log_line("PROGRESS task_id=", task_id, " scenario=", scenario_label,
             " replicate=", replicate_index, "/", task$reps,
             " seed=", seed, " state=START")
    out <- run_one_rw(seed = seed, scenario = scenario, n = task$n, m = task$m,
                      k = task$k, sim_root = sim_root)
    log_line("PROGRESS task_id=", task_id, " scenario=", scenario_label,
             " replicate=", replicate_index, "/", task$reps,
             " seed=", seed, " state=DONE")
    cbind(task[rep(1L, nrow(out)), c("task_id", "run_id", "chunk_id")],
          scenario_label = scenario_label, out)
  }))
}))
temporary_path <- tempfile(pattern = "rw_task_", tmpdir = raw_dir)
write.csv(raw, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, raw_path)) {
  unlink(temporary_path)
  stop("Could not publish task output.")
}
log_line("COMPLETE run_id=", task$run_id, " task_id=", task_id,
         " raw_path=", normalizePath(raw_path))
