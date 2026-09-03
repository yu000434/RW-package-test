#!/usr/bin/env Rscript

args <- commandArgs()
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
sim_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(sim_root, "R", "pmmrw.R"))
source(file.path(sim_root, "RW_scenario", "PMM", "scripts", "generate_rw_data.R"))

data <- make_rw_data(14000001L, n = 150L, type = "robins_2_2")$dat
method <- mice::make.method(data)
method[] <- ""
method["Z"] <- "pmm"
predictor_matrix <- mice::make.predictorMatrix(data)
predictor_matrix[,] <- 0L
predictor_matrix["Z", c("X", "A")] <- 1L

set.seed(14000002L)
ordinary <- mice::mice(
  data, m = 3L, method = method, predictorMatrix = predictor_matrix,
  blots = list(Z = list(donors = 5L)), print = FALSE
)
method["Z"] <- "pmmrw"
set.seed(14000002L)
recorded <- mice::mice(
  data, m = 3L, method = method, predictorMatrix = predictor_matrix,
  tasks = "train", blots = list(Z = list(donors = 5L)), print = FALSE
)

for (p in seq_len(ordinary$m)) {
  stopifnot(identical(mice::complete(ordinary, p)$Z, mice::complete(recorded, p)$Z))
}
donor_id <- rw::extract_donor_id(recorded, "Z")
missing <- which(is.na(data$Z))
for (p in seq_len(recorded$m)) {
  completed <- mice::complete(recorded, p)
  stopifnot(all(completed$Z[missing] == data$Z[donor_id[missing, p]]))
}

cat("MICE_COMPLETED_VALUE_PARITY=PASS\n")
cat("DONOR_ID_REPRODUCTION=PASS\n")
