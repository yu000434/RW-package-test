#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript audit_results.R TASK_FILE RAW_DIR [CELL_ID]")
}
tasks <- read.csv(args[[1L]], stringsAsFactors = FALSE)
raw_dir <- args[[2L]]
if (length(args) == 3L) {
  cell_id <- as.integer(args[[3L]])
  tasks <- tasks[tasks$cell_id == cell_id, , drop = FALSE]
  if (nrow(tasks) == 0L) stop("No tasks found for cell_id ", cell_id)
}
paths <- file.path(raw_dir, sprintf("gs_task%04d.csv", tasks$task_id))
if (any(!file.exists(paths))) stop("Missing task output: ", paste(basename(paths[!file.exists(paths)]), collapse = ", "))
raw <- do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
required <- c("task_id", "run_id", "cell_id", "chunk_id", "seed", "scenario", "n", "validated",
              "m", "k", "rb_total_var", "rr_total_var", "max_probability_error",
              "max_derivative_error", "max_intercept_error")
if (!all(required %in% names(raw))) stop("A required GS result column is missing.")
for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, , drop = FALSE]
  x <- raw[raw$task_id == task$task_id, , drop = FALSE]
  expected <- task$base_seed + seq_len(task$reps) - 1L
  if (nrow(x) != task$reps || !identical(sort(as.numeric(x$seed)), sort(as.numeric(expected)))) {
    stop("Task-level seed audit failed for task ", task$task_id)
  }
  if (any(x$run_id != task$run_id) || any(x$cell_id != task$cell_id) ||
      any(x$chunk_id != task$chunk_id) || any(x$scenario != "GS") ||
      any(x$n != task$n) || any(x$validated != task$validated) ||
      any(x$m != task$m) || any(x$k != task$k)) {
    stop("Task-level metadata audit failed for task ", task$task_id)
  }
}
expected_seed <- unlist(Map(function(base, reps) base + seq_len(reps) - 1L, tasks$base_seed, tasks$reps))
if (!identical(sort(as.numeric(raw$seed)), sort(as.numeric(expected_seed))) || anyDuplicated(raw$seed)) {
  stop("Seed audit failed.")
}
if (nrow(raw) != sum(tasks$reps)) stop("Replicate count audit failed.")
numeric_columns <- vapply(raw, is.numeric, logical(1))
required_finite <- setdiff(names(raw)[numeric_columns], c("validated", "quadrature_order"))
if (!all(is.finite(as.matrix(raw[, required_finite, drop = FALSE])))) stop("Non-finite numeric output.")
if (any(raw$rb_total_var <= 0) || any(raw$rr_total_var <= 0)) stop("Non-positive variance.")
if (max(raw$max_probability_error) > 1e-7 || max(raw$max_derivative_error) > 1e-7 ||
    max(raw$max_intercept_error) > 1e-7) stop("Soft-top-k identity audit failed.")
cat("RESULT_AUDIT=PASS\n")
cat("RUN_ID=", tasks$run_id[[1L]], " REPLICATES=", nrow(raw), "\n", sep = "")
