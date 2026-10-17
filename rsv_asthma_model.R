# Model functions. Source this file. It does nothing on its own.
# Everything is vectorized, so one number gives the base case and a
# vector of 10,000 gives the simulation.


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


# The same quantity as the first analysis computed it, in two passes.
# Kept so the first results can be reproduced.
prevalence_unexposed_published <- function(odds_ratio, prev_total,
                                           p_exposed) {
  start      <- prev_total / (odds_ratio * p_exposed + (1 - p_exposed))
  risk_ratio <- odds_ratio / ((1 - start) + start * odds_ratio)
  p0         <- prev_total / (risk_ratio * p_exposed + (1 - p_exposed))
  list(p0 = p0, p1 = p0 * risk_ratio, risk_ratio = risk_ratio)
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


# Median and 95% interval of a vector of trials.
interval <- function(x) {
  c(median = median(x),
    lower  = unname(quantile(x, 0.025)),
    upper  = unname(quantile(x, 0.975)))
}
