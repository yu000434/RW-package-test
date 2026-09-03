#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript make_draft_style_report.R RESULTS_DIR REPORTS_DIR")
}

results_dir <- normalizePath(args[[1L]])
reports_dir <- args[[2L]]
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

read_summary <- function(pattern) {
  paths <- list.files(results_dir, pattern = "^summary[.]csv$", full.names = TRUE, recursive = TRUE)
  paths <- paths[grepl(pattern, basename(dirname(paths)))]
  if (!length(paths)) stop("No matching summaries found for pattern: ", pattern)
  do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
}

gs <- read_summary("^gs")
gs <- gs[order(gs$m, gs$k), , drop = FALSE]
if (nrow(gs) != 16L || !all(sort(unique(gs$m)) == c(5L, 25L, 50L, 100L)) ||
    !all(sort(unique(gs$k)) == c(5L, 10L, 20L, 30L)) || any(gs$reps != 2500L)) {
  stop("The GS report requires the complete 4 by 4 grid with 2,500 replications per cell.")
}

rw_n1000 <- read.csv(file.path(results_dir, "rw1000_50_5", "summary.csv"), stringsAsFactors = FALSE)
rw_n150 <- read.csv(file.path(results_dir, "rw150_20_5", "summary.csv"), stringsAsFactors = FALSE)
scenario_order <- c("robins_1", "robins_2_1", "robins_2_2", "robins_2_3", "robins_3_1", "robins_3_2")
if (!identical(sort(rw_n1000$scenario), sort(scenario_order)) ||
    !identical(sort(rw_n150$scenario), sort(scenario_order)) ||
    any(rw_n1000$reps != 2500L) || any(rw_n150$reps != 2500L)) {
  stop("The RW report requires six scenarios with 2,500 replications each at both sample sizes.")
}
rw_n1000 <- rw_n1000[match(scenario_order, rw_n1000$scenario), , drop = FALSE]
rw_n150 <- rw_n150[match(scenario_order, rw_n150$scenario), , drop = FALSE]

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required.")
long <- rbind(
  data.frame(m = gs$m, k = factor(gs$k, levels = c(5, 10, 20, 30)),
             coverage_type = "A) Empirical-centered coverage", method = "Rubin-PMM",
             coverage = gs$rr_coverage_empirical_center,
             lower = gs$rr_coverage_empirical_center_lower,
             upper = gs$rr_coverage_empirical_center_upper),
  data.frame(m = gs$m, k = factor(gs$k, levels = c(5, 10, 20, 30)),
             coverage_type = "A) Empirical-centered coverage", method = "PMM-corrected",
             coverage = gs$rb_coverage_empirical_center,
             lower = gs$rb_coverage_empirical_center_lower,
             upper = gs$rb_coverage_empirical_center_upper),
  data.frame(m = gs$m, k = factor(gs$k, levels = c(5, 10, 20, 30)),
             coverage_type = "B) True-parameter coverage", method = "Rubin-PMM",
             coverage = gs$rr_coverage_beta0,
             lower = gs$rr_coverage_beta0_lower,
             upper = gs$rr_coverage_beta0_upper),
  data.frame(m = gs$m, k = factor(gs$k, levels = c(5, 10, 20, 30)),
             coverage_type = "B) True-parameter coverage", method = "PMM-corrected",
             coverage = gs$rb_coverage_beta0,
             lower = gs$rb_coverage_beta0_lower,
             upper = gs$rb_coverage_beta0_upper)
)
long$method <- factor(long$method, levels = c("Rubin-PMM", "PMM-corrected"))
long$coverage_type <- factor(long$coverage_type,
                             levels = c("A) Empirical-centered coverage", "B) True-parameter coverage"))
figure_path <- file.path(reports_dir, "coverage.pdf")
plot <- ggplot2::ggplot(long, ggplot2::aes(x = m, y = coverage, group = method,
                                           linetype = method, shape = method)) +
  ggplot2::geom_hline(yintercept = 0.95, colour = "grey40", linewidth = 0.35) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 1.8, linewidth = 0.35) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::geom_point(size = 2.1, fill = "white") +
  ggplot2::facet_grid(
    coverage_type ~ k,
    labeller = ggplot2::labeller(
      k = function(x) paste("k =", x),
      coverage_type = function(x) x
    )
  ) +
  ggplot2::scale_x_continuous(breaks = c(5, 25, 50, 100)) +
  ggplot2::scale_y_continuous(breaks = c(0.85, 0.90, 0.95, 1.00), limits = c(0.83, 1.01)) +
  ggplot2::scale_linetype_manual(values = c("Rubin-PMM" = "dashed", "PMM-corrected" = "solid")) +
  ggplot2::scale_shape_manual(values = c("Rubin-PMM" = 0, "PMM-corrected" = 16)) +
  ggplot2::labs(x = "Number of imputations (m)", y = "Coverage", linetype = NULL, shape = NULL) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey45"),
    panel.spacing = grid::unit(0.13, "in"),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
  )
ggplot2::ggsave(figure_path, plot, width = 10.2, height = 5.8)

num <- function(x) formatC(x, format = "f", digits = 3L)
pair <- function(x, y) paste(num(x), num(y), sep = " / ")
gs_rows <- unlist(lapply(unique(gs$m), function(m) {
  x <- gs[gs$m == m, , drop = FALSE]
  rows <- vapply(seq_len(nrow(x)), function(i) {
    paste(
      if (i == 1L) x$m[i] else "",
      sprintf("PMM ($k=%d$)", x$k[i]),
      num(x$mean_estimate[i]), num(x$bias[i]), num(x$empirical_sd[i]),
      pair(x$mean_rb_se[i], x$mean_rr_se[i]),
      pair(x$rb_coverage_empirical_center[i], x$rr_coverage_empirical_center[i]),
      pair(x$rb_coverage_beta0[i], x$rr_coverage_beta0[i]),
      sep = " & "
    ) |> paste0(" \\\\")
  }, character(1))
  if (m != tail(unique(gs$m), 1L)) rows <- c(rows, "\\addlinespace[2pt]")
  rows
}), use.names = FALSE)

scenario_labels <- c(robins_1 = "R1", robins_2_1 = "R2a", robins_2_2 = "R2b",
                     robins_2_3 = "R2c", robins_3_1 = "R3a", robins_3_2 = "R3b")
rw_rows <- function(x) vapply(seq_len(nrow(x)), function(i) {
  paste(
    scenario_labels[[x$scenario[i]]], "PMM ($k=5$)",
    num(x$mean_estimate[i]), num(x$bias[i]), num(x$empirical_sd[i]),
    pair(x$mean_rb_se[i], x$mean_rr_se[i]),
    pair(x$rb_coverage_empirical_center[i], x$rr_coverage_empirical_center[i]),
    pair(x$rb_coverage_beta0[i], x$rr_coverage_beta0[i]),
    sep = " & "
  ) |> paste0(" \\\\")
}, character(1))

table_header_gs <- c(
  "\\begin{tabular}{rlrrrrrr}", "\\toprule",
  "$m$ & \\makecell{Imputation\\\\model} & \\makecell{Mean\\\\estimate} & Bias &",
  "\\makecell{Empirical\\\\SD} & \\makecell{Mean estimated SE\\\\Proposed/Rubin} &",
  "\\makecell{Coverage $\\beta^\\ast$\\\\Proposed/Rubin} &",
  "\\makecell{Coverage $\\beta_0$\\\\Proposed/Rubin} \\\\", "\\midrule"
)
table_header_rw <- c(
  "\\begin{tabular}{llrrrrrr}", "\\toprule",
  "Scenario & \\makecell{Imputation\\\\model} & \\makecell{Mean\\\\estimate} & Bias &",
  "\\makecell{Empirical\\\\SD} & \\makecell{Mean estimated SE\\\\Proposed/Rubin} &",
  "\\makecell{Coverage $\\beta^\\ast$\\\\Proposed/Rubin} &",
  "\\makecell{Coverage $\\beta_0$\\\\Proposed/Rubin} \\\\", "\\midrule"
)

method_pages <- c(
  "\\begin{center}",
  "{\\Large MICE PMM Variance Estimation}\\\\[4pt]",
  "{\\normalsize Standard MICE PMM with a soft top-k donor-source correction}",
  "\\end{center}",
  "\\section*{1. Goal and notation}",
  "For each simulation replicate, we create $m$ completed data sets with ordinary \\texttt{mice::mice(method = \"pmm\")}. The completed values and realized donor IDs are exactly those produced by standard MICE PMM. If donor $j$ is selected for recipient $i$, the completed value for $i$ is the observed value of $j$.",
  "",
  "The variance problem is that one donor can be used for several recipients. The copied values then share the same observed donor source. A rowwise variance calculation misses this dependence. The correction therefore keeps the realized donor map when forming the completed-data score covariance.",
  "",
  "Let $U_i^{(p)}$ be the completed-data analysis score for person $i$ in imputation $p$, let $d_i^{(p)}$ be the imputation-model influence quantity, and let $\\tau^{(p)}$ be the completed-data analysis information. For each imputation, the realized donor IDs define a donor-source version of the score, $U_{\\Omega,i}^{(p)}$: every recipient who received donor $j$ is assigned to donor source $j$ when the score covariance is formed.",
  "",
  "The key donor-reuse component is",
  "\\[",
  "  \\widehat{\\Omega} = \\frac{1}{n}\\bar U_{\\Omega}^{\\mathsf T}\\bar U_{\\Omega},",
  "  \\qquad",
  "  \\bar U_{\\Omega} = \\frac{1}{m}\\sum_{p=1}^m U_{\\Omega}^{(p)}.",
  "\\]",
  "Thus, donor reuse affects $\\widehat{\\Omega}$ through the actual MICE donor labels, not through a simulated approximation to the completed values.",
  "",
  "\\section*{2. RW variance structure}",
  "The estimator follows the Robins--Wang form. With $\\bar d=m^{-1}\\sum_p d^{(p)}$, $\\widehat{\\alpha}=(nm)^{-1}\\sum_p d^{(p)\\mathsf T}d^{(p)}$, and $\\widehat{\\tau}=(nm)^{-1}\\sum_p\\tau^{(p)}$, define",
  "\\[",
  "  \\widehat{\\Delta} = \\widehat{\\Omega} + \\widehat{\\kappa}\\widehat{\\alpha}\\widehat{\\kappa}^{\\mathsf T}",
  "  + \\frac{\\widehat{\\kappa}\\bar d^{\\mathsf T}\\bar U_{\\Omega} + (\\widehat{\\kappa}\\bar d^{\\mathsf T}\\bar U_{\\Omega})^{\\mathsf T}}{n},",
  "  \\qquad",
  "  \\widehat V = \\frac{1}{n}\\widehat{\\tau}^{-1}\\widehat{\\Delta}\\widehat{\\tau}^{-\\mathsf T}.",
  "\\]",
  "The remaining task is to construct the PMM part of $\\widehat{\\kappa}$ in a way that reflects the matching rule while preserving standard MICE PMM draws.",
  "\\clearpage",
  "\\section*{3. Soft top-k working probability}",
  "For a missing recipient $i$ and an observed donor $j$, standard PMM selects uniformly from the recipient's nearest-$k$ donor set. That selection rule is discontinuous as the PMM regression coefficient $\\psi$ changes. We use a smooth probability only for the variance score.",
  "",
  "Let $l_j$ and $u_j$ be the lower and upper prediction-score boundaries for donor $j$ to enter a nearest-$k$ set. Let $\\eta_i$ be the recipient's PMM prediction and $s_i$ be the fitted posterior prediction standard deviation. The working probability is",
  "\\[",
  " q_{ij}(\\psi) = \\frac{1}{k}\\left[",
  " \\Phi\\left\\{\\frac{u_j-\\eta_i}{s_i}\\right\\} -",
  " \\Phi\\left\\{\\frac{l_j-\\eta_i}{s_i}\\right\\}",
  " \\right].",
  "\\]",
  "Its derivative moves the donor predictions, recipient prediction, and the two membership boundaries together. The probability is used to calculate a differentiable PMM score. It does not replace MICE matching: MICE still chooses one donor uniformly from its realized nearest-$k$ set and copies that donor's observed value.",
  "",
  "\\section*{4. Rao--Blackwell donor-source term}",
  "Let $U_{ij}^{(p)}$ denote the counterfactual analysis score if recipient $i$ were assigned donor $j$ in imputation $p$. Directly using the observed donor outcome in this term is noisy. We instead use its conditional expectation,",
  "\\[",
  " G_{ij}^{(p)} = E\\left\\{U_{ij}^{(p)} \\mid X_j, \\text{recipient information}\\right\\},",
  " \\qquad",
  " \\widehat{\\kappa}_{\\mathrm{RB}} = \\frac{1}{nm}\\sum_{p=1}^m\\sum_i\\sum_j",
  " G_{ij}^{(p)}\\frac{\\partial q_{ij}^{(p)}}{\\partial\\psi}.",
  "\\]",
  "For the RW linear setting, $G_{ij}^{(p)}$ is obtained by replacing the donor outcome with its fitted normal-model mean. For the GS setting, $G_{ij}^{(p)}$ integrates the downstream score under that fitted normal PMM model using fixed Gauss--Legendre quadrature.",
  "",
  "The derivative applies only to $q_{ij}$. We do not add $q_{ij}\\,\\partial G_{ij}/\\partial\\psi$, because MICE matching conditions on observed donor values; the variance score describes how donor selection changes, not a newly generated donor outcome.",
  "\\clearpage",
  "\\section*{5. Calculation sequence and reported comparisons}",
  "For each imputed data set, the calculation is:",
  "\\begin{enumerate}",
  "\\item Run standard MICE PMM and record its completed values and realized donor IDs.",
  "\\item Fit the completed-data analysis model and obtain $U^{(p)}$, $d^{(p)}$, and $\\tau^{(p)}$.",
  "\\item Construct $q_{ij}^{(p)}$ and $\\partial q_{ij}^{(p)}/\\partial\\psi$ from the fitted PMM prediction model.",
  "\\item Form $\\widehat{\\kappa}_{\\mathrm{RB}}$, retain the actual donor-source score in $\\widehat{\\Omega}$, and assemble $\\widehat V$ with the RW formula.",
  "\\end{enumerate}",
  "The tables report PMM-corrected/Rubin-PMM pairs. The point estimate is the same under both methods because both use the same completed MICE PMM data. They differ only in the estimated variance and hence in coverage. Coverage of $\\beta^\\ast$ uses the Monte Carlo mean estimate as the target and assesses variance calibration. Coverage of $\\beta_0$ uses the data-generating coefficient and therefore also reflects finite-sample bias.",
  "",
  "All completed cells below use 2,500 Monte Carlo replications. The figure shows exact 95\\% binomial Monte Carlo intervals for each coverage estimate.",
  "\\clearpage"
)

tex <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=0.7in]{geometry}",
  "\\usepackage{amsmath,amssymb,booktabs,makecell,graphicx,float}",
  "\\usepackage[T1]{fontenc}",
  "\\begin{document}",
  method_pages,
  "\\begin{center}",
  "{\\Large MICE PMM Simulation Results}\\\\[4pt]",
  "{\\normalsize Current standard-MICE PMM results; 2,500 Monte Carlo replications per completed cell}",
  "\\end{center}",
  "\\vspace{4pt}",
  "\\begin{table}[H]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{2.6pt}",
  table_header_gs, gs_rows, "\\bottomrule", "\\end{tabular}",
  "\\caption{GS setting, $n=2000$. Proposed/Rubin pairs are PMM-corrected/Rubin-PMM. Bias is the Monte Carlo mean estimate minus the data-generating coefficient.}",
  "\\end{table}", "\\clearpage",
  "\\begin{figure}[H]", "\\centering",
  "\\includegraphics[width=\\textwidth]{coverage.pdf}",
  "\\caption{GS coverage for standard MICE PMM. Error bars are exact 95\\% binomial Monte Carlo intervals.}",
  "\\end{figure}", "\\clearpage",
  "\\begin{table}[H]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3.2pt}",
  table_header_rw, rw_rows(rw_n1000), "\\bottomrule", "\\end{tabular}",
  "\\caption{RW setting, $n=1000$, $m=50$, and $k=5$. Each scenario has 2,500 Monte Carlo replications. Proposed/Rubin pairs are PMM-corrected/Rubin-PMM.}",
  "\\end{table}", "\\clearpage",
  "\\begin{table}[H]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3.2pt}",
  table_header_rw, rw_rows(rw_n150), "\\bottomrule", "\\end{tabular}",
  "\\caption{RW setting, $n=150$, $m=20$, and $k=5$. Each scenario has 2,500 Monte Carlo replications. Proposed/Rubin pairs are PMM-corrected/Rubin-PMM.}",
  "\\end{table}", "\\end{document}"
)
tex_path <- file.path(reports_dir, "results.tex")
writeLines(tex, tex_path)
cat("FIGURE_PDF=", normalizePath(figure_path), "\n", sep = "")
cat("REPORT_TEX=", normalizePath(tex_path), "\n", sep = "")
