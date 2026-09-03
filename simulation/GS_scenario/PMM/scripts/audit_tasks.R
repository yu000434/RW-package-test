#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript audit_tasks.R TASK_FILE")
task_file <- args[[1]]
tasks <- read.csv(task_file, stringsAsFactors = FALSE)
required <- c("task_id", "run_id", "cell_id", "scenario", "n", "validated", "m", "k",
              "chunk_id", "chunks_per_cell", "reps", "reps_per_cell", "base_seed")
if (!identical(names(tasks), required)) stop("Task table columns do not match the GS schema.")
if (!identical(tasks$task_id, seq_len(nrow(tasks)))) stop("Task IDs are not consecutive.")
if (length(unique(tasks$run_id)) != 1L || any(tasks$scenario != "GS")) stop("Run ID or scenario is invalid.")

cell_key <- interaction(tasks$cell_id, tasks$m, tasks$k, drop = TRUE)
for (key in levels(cell_key)) {
  x <- tasks[cell_key == key, , drop = FALSE]
  x <- x[order(x$chunk_id), , drop = FALSE]
  if (!identical(x$chunk_id, seq_len(nrow(x)))) stop("Chunk IDs are invalid for cell ", key)
  if (any(x$chunks_per_cell != nrow(x))) stop("Chunk count is invalid for cell ", key)
  if (sum(x$reps) != x$reps_per_cell[1L]) stop("Replicate total is invalid for cell ", key)
  starts <- x$base_seed
  expected <- starts[1L] + cumsum(c(0L, head(x$reps, -1L)))
  if (!identical(starts, expected)) stop("Seed blocks are not consecutive for cell ", key)
}
all_seeds <- unlist(Map(function(base, reps) base + seq_len(reps) - 1L, tasks$base_seed, tasks$reps))
if (anyDuplicated(all_seeds)) stop("Seeds overlap across GS cells.")
cat("TASK_AUDIT=PASS\n")
cat("RUN_ID=", tasks$run_id[1L], " CELLS=", length(levels(cell_key)),
    " TASKS=", nrow(tasks), " REPLICATES=", sum(tasks$reps), "\n", sep = "")
