#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- sub("^--[^=]+=", "", args)
  names(values) <- sub("^--([^=]+)=.*$", "\\1", args)
  values
}

scenarios <- c(S1 = "robins_1", S2a = "robins_2_1", S2b = "robins_2_2",
               S2c = "robins_2_3", S3a = "robins_3_1", S3b = "robins_3_2")
seed_block <- 100000L
required <- c("n", "m", "k", "reps", "chunk-reps", "seed-base", "run-id", "task-file")
args <- parse_args(commandArgs(trailingOnly = TRUE))
missing <- setdiff(required, names(args))
if (length(missing)) {
  stop("Missing arguments: ", paste(paste0("--", missing, "="), collapse = ", "))
}

n <- as.integer(args[["n"]])
m <- as.integer(args[["m"]])
k <- as.integer(args[["k"]])
reps_per_scenario <- as.integer(args[["reps"]])
chunk_reps <- as.integer(args[["chunk-reps"]])
seed_base <- as.integer(args[["seed-base"]])
run_id <- args[["run-id"]]
task_file <- args[["task-file"]]
if (anyNA(c(n, m, k, reps_per_scenario, chunk_reps, seed_base)) ||
    n <= 0L || m <= 0L || k <= 0L || reps_per_scenario <= 0L ||
    chunk_reps <= 0L || seed_base < 1L) {
  stop("Numeric arguments are invalid.")
}
if (reps_per_scenario > seed_block) {
  stop("`reps` must not exceed 100000 because each cell has a 100000-seed block.")
}
if (!nzchar(run_id)) stop("`run-id` must not be empty.")
if (!grepl("^[A-Za-z0-9_-]+$", run_id)) {
  stop("`run-id` may contain only letters, digits, underscores, and hyphens.")
}
if (file.exists(task_file)) stop("Task file already exists: ", task_file)

max_seed <- seed_base + (length(scenarios) - 1L) * seed_block + reps_per_scenario - 1L
if (max_seed > .Machine$integer.max) stop("The requested seed range exceeds R's integer range.")
n_chunks <- ceiling(reps_per_scenario / chunk_reps)
reps <- rep.int(chunk_reps, n_chunks)
reps[n_chunks] <- reps_per_scenario - chunk_reps * (n_chunks - 1L)
tasks <- data.frame(
  task_id = seq_len(n_chunks), run_id = run_id, n = n, m = m, k = k,
  chunk_id = seq_len(n_chunks), chunks_per_run = n_chunks,
  reps = reps, reps_per_scenario = reps_per_scenario,
  base_seed = seed_base + cumsum(c(1L, head(reps, -1L))) - 1L
)

dir.create(dirname(task_file), recursive = TRUE, showWarnings = FALSE)
write.csv(tasks, task_file, row.names = FALSE)
cat("WROTE_TASK_FILE=", normalizePath(task_file, mustWork = TRUE), "\n", sep = "")
cat("SCENARIOS=", length(scenarios), " TASKS=", nrow(tasks),
    " REPLICATES_PER_SCENARIO=", sum(tasks$reps), "\n", sep = "")
