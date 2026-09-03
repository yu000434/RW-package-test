#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- sub("^--[^=]+=", "", args)
  names(values) <- sub("^--([^=]+)=.*$", "\\1", args)
  values
}

int_values <- function(value, name) {
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (!length(out) || anyNA(out) || any(out <= 0L)) {
    stop("`", name, "` must be a comma-separated list of positive integers.")
  }
  out
}

required <- c("n", "validated", "m", "k", "reps", "chunk-reps", "seed-base", "run-id", "task-file")
args <- parse_args(commandArgs(trailingOnly = TRUE))
missing <- setdiff(required, names(args))
if (length(missing)) {
  stop("Missing arguments: ", paste(paste0("--", missing, "="), collapse = ", "))
}

n <- as.integer(args[["n"]])
validated <- as.integer(args[["validated"]])
m_values <- int_values(args[["m"]], "m")
k_values <- int_values(args[["k"]], "k")
reps_per_cell <- as.integer(args[["reps"]])
chunk_reps <- as.integer(args[["chunk-reps"]])
seed_base <- as.integer(args[["seed-base"]])
run_id <- args[["run-id"]]
task_file <- args[["task-file"]]

if (anyNA(c(n, validated, reps_per_cell, chunk_reps, seed_base)) ||
    n <= 0L || validated <= 0L || validated >= n || reps_per_cell <= 0L ||
    chunk_reps <= 0L || seed_base < 1L) {
  stop("`n`, `validated`, `reps`, `chunk-reps`, and `seed-base` are invalid.")
}
if (reps_per_cell > 100000L) {
  stop("`reps` must not exceed 100000 because each cell has a 100000-seed block.")
}
if (!nzchar(run_id)) stop("`run-id` must not be empty.")
if (!grepl("^[A-Za-z0-9_-]+$", run_id)) {
  stop("`run-id` may contain only letters, digits, underscores, and hyphens.")
}
if (file.exists(task_file)) stop("Task file already exists: ", task_file)

cells <- expand.grid(m = m_values, k = k_values, KEEP.OUT.ATTRS = FALSE)
cells <- cells[order(cells$m, cells$k), , drop = FALSE]
if (anyDuplicated(cells[c("m", "k")])) stop("Each GS (m, k) cell must occur once.")
cells$cell_id <- seq_len(nrow(cells))
max_seed <- seed_base + (nrow(cells) - 1L) * 100000L + reps_per_cell - 1L
if (max_seed > .Machine$integer.max) stop("The requested seed range exceeds R's integer range.")

tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  n_chunks <- ceiling(reps_per_cell / chunk_reps)
  reps <- rep.int(chunk_reps, n_chunks)
  reps[n_chunks] <- reps_per_cell - chunk_reps * (n_chunks - 1L)
  data.frame(
    run_id = run_id,
    cell_id = cells$cell_id[i],
    scenario = "GS",
    n = n,
    validated = validated,
    m = cells$m[i],
    k = cells$k[i],
    chunk_id = seq_len(n_chunks),
    chunks_per_cell = n_chunks,
    reps = reps,
    reps_per_cell = reps_per_cell,
    base_seed = seed_base + (cells$cell_id[i] - 1L) * 100000L + cumsum(c(1L, head(reps, -1L))) - 1L
  )
}))
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[c("task_id", "run_id", "cell_id", "scenario", "n", "validated", "m", "k",
                 "chunk_id", "chunks_per_cell", "reps", "reps_per_cell", "base_seed")]

dir.create(dirname(task_file), recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, task_file, row.names = FALSE)
cat("WROTE_TASK_FILE=", normalizePath(task_file, mustWork = TRUE), "\n", sep = "")
cat("CELLS=", nrow(cells), " TASKS=", nrow(tasks), " REPLICATES=", sum(tasks$reps), "\n", sep = "")
