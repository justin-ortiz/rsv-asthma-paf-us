# compare_calibrations.R
#
# Prints the base case both ways: as the original analysis computed
# prevalence in unexposed children, and as the Methods describes it.
#
# Run this file on its own. It changes nothing.

source("rsv_asthma_model.R")

pop_total  <- 3664292
prev_total <- 0.088
odds_ratio <- 3.02

inputs <- readRDS("inputs.rds")

correction <- 0.491 / (0.491 + 0.374) * 69068  / 107338 +
              0.374 / (0.491 + 0.374) * 170744 / 256495

n_natural    <- unname(inputs$scenarios[["Natural History"]]$point) * correction
n_nirsevimab <- unname(inputs$scenarios[["Nirsevimab"]]$point)      * correction

for (method in c("published", "closed form")) {

  p <- if (method == "published")
         prevalence_unexposed_published(odds_ratio, prev_total,
                                        n_natural / pop_total)
       else
         prevalence_unexposed(odds_ratio, prev_total,
                              n_natural / pop_total)

  natural    <- outcomes(n_natural,    pop_total, p$p0, p$risk_ratio)
  nirsevimab <- outcomes(n_nirsevimab, pop_total, p$p0, p$risk_ratio)
  prevented  <- natural$attributable - nirsevimab$attributable

  cat("\n", method, "\n", sep = "")
  cat("  prevalence, unexposed      ", round(p$p0, 5), "\n")
  cat("  prevalence, exposed        ", round(p$p1, 5), "\n")
  cat("  relative risk              ", round(p$risk_ratio, 3), "\n")
  cat("  all-cause asthma cases     ", round(natural$asthma_total), "\n")
  cat("  cases if none exposed      ", round(natural$asthma_if_none_exposed), "\n")
  cat("  attributable cases         ", round(natural$attributable), "\n")
  cat("  attributable fraction      ",
      round(100 * natural$attributable_fraction, 2), "%\n", sep = "")
  cat("  prevented by nirsevimab    ", round(prevented), "\n")
  cat("  prevented, % of attributable",
      round(100 * prevented / natural$attributable, 1), "\n")
}
