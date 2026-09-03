#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript audit_tasks.R TASK_FILE")
task_file <- args[[1]]
tasks <- read.csv(task_file, stringsAsFactors = FALSE)
scenario_count <- 6L
required <- c("task_id", "run_id", "n", "m", "k", "chunk_id", "chunks_per_run",
              "reps", "reps_per_scenario", "base_seed")
if (!identical(names(tasks), required)) stop("Task table columns do not match the RW schema.")
if (!identical(tasks$task_id, seq_len(nrow(tasks)))) stop("Task IDs are not consecutive.")
if (length(unique(tasks$run_id)) != 1L) stop("Run ID is invalid.")

tasks <- tasks[order(tasks$chunk_id), , drop = FALSE]
if (!identical(tasks$chunk_id, seq_len(nrow(tasks)))) stop("Chunk IDs are invalid.")
if (any(tasks$chunks_per_run != nrow(tasks))) stop("Chunk count is invalid.")
if (sum(tasks$reps) != tasks$reps_per_scenario[1L]) stop("Replicate total is invalid.")
expected <- tasks$base_seed[1L] + cumsum(c(0L, head(tasks$reps, -1L)))
if (!identical(tasks$base_seed, expected)) stop("Seed blocks are not consecutive.")
all_seeds <- unlist(Map(function(base, reps) base + seq_len(reps) - 1L, tasks$base_seed, tasks$reps))
all_scenario_seeds <- unlist(lapply(0:(scenario_count - 1L), function(index) {
  all_seeds + index * 100000L
}))
if (anyDuplicated(all_scenario_seeds)) stop("Seeds overlap across scenarios.")
cat("TASK_AUDIT=PASS\n")
cat("RUN_ID=", tasks$run_id[1L], " SCENARIOS=", scenario_count,
    " TASKS=", nrow(tasks), " REPLICATES_PER_SCENARIO=", sum(tasks$reps), "\n", sep = "")
