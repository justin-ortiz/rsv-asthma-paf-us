# Functions used by run_analysis.R.


# Asthma prevalence at age 6 in children without infant RSV LRTI,
# solved from national prevalence and the odds ratio. Closed form, as
# in the Methods.
prevalence_unexposed <- function(odds_ratio, prev_total, p_exposed) {
  a  <- (1 - p_exposed) * (odds_ratio - 1)
  b  <- p_exposed * odds_ratio + (1 - p_exposed) -
        prev_total * (odds_ratio - 1)
  p0 <- (-b + sqrt(b^2 + 4 * a * prev_total)) / (2 * a)
  p1 <- odds_ratio * p0 / (1 - p0 + odds_ratio * p0)
  list(p0 = p0, p1 = p1, risk_ratio = p1 / p0)
}


# Asthma outcomes for one scenario. Calibrate once on natural history
# and pass the same p0 and risk ratio to every scenario.
outcomes <- function(n_exposed, pop_total, p0, risk_ratio) {
  asthma_total <- (pop_total - n_exposed) * p0 + n_exposed * p0 * risk_ratio
  attributable <- n_exposed * p0 * (risk_ratio - 1)
  data.frame(
    n_exposed              = n_exposed,
    p_exposed              = n_exposed / pop_total,
    prev_unexposed         = p0,
    prev_exposed           = p0 * risk_ratio,
    risk_ratio             = risk_ratio,
    asthma_total           = asthma_total,
    asthma_if_none_exposed = pop_total * p0,
    attributable           = attributable,
    attributable_fraction  = attributable / asthma_total)
}


# 95% uncertainty interval.
interval <- function(x) {
  c(lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975)))
}
