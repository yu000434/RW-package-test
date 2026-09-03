#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript audit_results.R TASK_FILE RAW_DIR")
}
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
raw_dir <- args[[2L]]
scenarios <- c(S1 = "robins_1", S2a = "robins_2_1", S2b = "robins_2_2",
               S2c = "robins_2_3", S3a = "robins_3_1", S3b = "robins_3_2")
seed_block <- 100000L
paths <- file.path(raw_dir, sprintf("rw_task%04d.csv", tasks$task_id))
if (any(!file.exists(paths))) stop("Missing task output: ", paste(basename(paths[!file.exists(paths)]), collapse = ", "))
raw <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
required <- c("task_id", "run_id", "chunk_id", "scenario_label", "seed", "scenario",
              "n", "m", "k", "rb_total_var", "rr_total_var", "max_probability_error",
              "max_derivative_error", "max_intercept_error")
if (!all(required %in% names(raw))) stop("A required RW result column is missing.")
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  x <- raw[raw$task_id == task$task_id, , drop = FALSE]
  if (nrow(x) != task$reps * length(scenarios)) stop("Task row count audit failed for task ", task$task_id)
  if (any(x$run_id != task$run_id) || any(x$chunk_id != task$chunk_id) || any(x$n != task$n) ||
      any(x$m != task$m) || any(x$k != task$k)) {
    stop("Task-level metadata audit failed for task ", task$task_id)
  }
  for (scenario_index in seq_along(scenarios)) {
    scenario_label <- names(scenarios)[[scenario_index]]
    y <- x[x$scenario_label == scenario_label, , drop = FALSE]
    expected <- task$base_seed + seq_len(task$reps) - 1L + (scenario_index - 1L) * seed_block
    if (nrow(y) != task$reps || !identical(sort(as.numeric(y$seed)), sort(as.numeric(expected))) ||
        any(y$scenario != scenarios[[scenario_index]])) {
      stop("Scenario seed audit failed for task ", task$task_id, ", scenario ", scenario_label)
    }
  }
}
base_seed <- unlist(Map(function(base, reps) base + seq_len(reps) - 1L, tasks$base_seed, tasks$reps))
expected_seed <- unlist(lapply(0:(length(scenarios) - 1L), function(index) base_seed + index * seed_block))
if (!identical(sort(as.numeric(raw$seed)), sort(as.numeric(expected_seed))) || anyDuplicated(raw$seed)) {
  stop("Seed audit failed.")
}
if (nrow(raw) != sum(tasks$reps) * length(scenarios)) stop("Replicate count audit failed.")
numeric_columns <- vapply(raw, is.numeric, logical(1))
required_finite <- setdiff(names(raw)[numeric_columns], c("validated", "quadrature_order"))
if (!all(is.finite(as.matrix(raw[, required_finite, drop = FALSE])))) stop("Non-finite numeric output.")
if (any(raw$rb_total_var <= 0) || any(raw$rr_total_var <= 0)) stop("Non-positive variance.")
if (max(raw$max_probability_error) > 1e-7 || max(raw$max_derivative_error) > 1e-7 ||
    max(raw$max_intercept_error) > 1e-7) stop("Soft-top-k identity audit failed.")
cat("RESULT_AUDIT=PASS\n")
cat("RUN_ID=", tasks$run_id[[1L]], " SCENARIOS=", length(scenarios),
    " REPLICATES_PER_SCENARIO=", sum(tasks$reps), "\n", sep = "")
