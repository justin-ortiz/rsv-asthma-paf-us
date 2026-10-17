# Edit the assumptions, then run the whole file.

source("rsv_asthma_model.R")


# Assumptions

set.seed(20261017)
trials <- 10000

input_cache <- "inputs.rds"
scenario    <- "Natural History"

pop_total <- 3664292          # NCHS 2021 births

prev_total       <- 0.088     # NHIS 2018-2021, asthma at age 6
prev_total_lower <- 0.081
prev_total_upper <- 0.096

episode_private  <- 69068     # Gantenberg, visit to episode correction
visit_private    <- 107338
episode_medicaid <- 170744
visit_medicaid   <- 256495
share_private    <- 0.491
share_medicaid   <- 0.374

# Row 1 is the base case. Add a row to add an analysis.
effects <- data.frame(
  label      = c("primary",
                 "published primary model",
                 "unconditional, asthma at age 6 or later",
                 "full dataset with three indicators",
                 "unconditional, all three conditions met"),
  odds_ratio = c(3.02, 2.45, 2.64, 2.91, 2.95),
  lower      = c(1.79, 1.23, 1.75, 1.77, 1.65),
  upper      = c(5.11, 4.88, 3.98, 4.78, 5.28))


# Inputs

inputs     <- readRDS(input_cache)
encounters <- inputs$scenarios[[scenario]]

w_private  <- share_private  / (share_private + share_medicaid)
w_medicaid <- share_medicaid / (share_private + share_medicaid)

correction_point <- w_private  * episode_private  / visit_private +
                    w_medicaid * episode_medicaid / visit_medicaid

correction_trials <-
  w_private  * rbeta(trials, episode_private + 1,
                     visit_private - episode_private + 1) +
  w_medicaid * rbeta(trials, episode_medicaid + 1,
                     visit_medicaid - episode_medicaid + 1)

n_exposed_point  <- encounters$point  * correction_point
n_exposed_trials <- encounters$trials * correction_trials

prev_total_trials <- rnorm(trials, prev_total,
                           (prev_total_upper - prev_total_lower) / (2 * 1.96))

# One set of normal draws, reused for every odds ratio, so rows differ
# only by the effect estimate.
z <- rnorm(trials)


# Every effect estimate, base case and simulation

results <- data.frame()

for (i in 1:nrow(effects)) {

  log_sd     <- (log(effects$upper[i]) - log(effects$lower[i])) / (2 * 1.96)
  odds_draws <- exp(log(effects$odds_ratio[i]) + log_sd * z)

  cal_point <- prevalence_unexposed(effects$odds_ratio[i], prev_total,
                                    n_exposed_point / pop_total)
  base <- outcomes(n_exposed_point, pop_total,
                   cal_point$p0, cal_point$risk_ratio)

  cal_sim <- prevalence_unexposed(odds_draws, prev_total_trials,
                                  n_exposed_trials / pop_total)
  sim <- outcomes(n_exposed_trials, pop_total,
                  cal_sim$p0, cal_sim$risk_ratio)

  fraction <- interval(sim$attributable_fraction)
  cases    <- interval(sim$attributable)

  results <- rbind(results, data.frame(
    label          = effects$label[i],
    odds_ratio     = effects$odds_ratio[i],
    risk_ratio     = round(base$risk_ratio, 3),
    prev_unexposed = round(base$prev_unexposed, 5),
    fraction       = round(100 * base$attributable_fraction, 2),
    fraction_lower = round(100 * fraction["lower"], 1),
    fraction_upper = round(100 * fraction["upper"], 1),
    cases          = round(base$attributable),
    cases_lower    = round(cases["lower"]),
    cases_upper    = round(cases["upper"])))
}

cat("\nEFFECT ESTIMATES\n\n")
print(results, row.names = FALSE)


# Scenarios against natural history, base case odds ratio.
# Calibrated once on natural history and held fixed. Differences are
# taken trial by trial, which works because all scenarios in one
# workbook share the same trials.

log_sd     <- (log(effects$upper[1]) - log(effects$lower[1])) / (2 * 1.96)
odds_draws <- exp(log(effects$odds_ratio[1]) + log_sd * z)

cal_point <- prevalence_unexposed(effects$odds_ratio[1], prev_total,
                                  n_exposed_point / pop_total)
cal_sim   <- prevalence_unexposed(odds_draws, prev_total_trials,
                                  n_exposed_trials / pop_total)

nh_point <- outcomes(n_exposed_point,  pop_total,
                     cal_point$p0, cal_point$risk_ratio)
nh_sim   <- outcomes(n_exposed_trials, pop_total,
                     cal_sim$p0, cal_sim$risk_ratio)

scenarios <- data.frame()

for (s in c("Nirsevimab", "RSVpreF")) {

  s_point <- outcomes(unname(inputs$scenarios[[s]]$point) * correction_point,
                      pop_total, cal_point$p0, cal_point$risk_ratio)
  s_sim   <- outcomes(unname(inputs$scenarios[[s]]$trials) * correction_trials,
                      pop_total, cal_sim$p0, cal_sim$risk_ratio)

  prevented <- nh_point$attributable - s_point$attributable
  spread    <- interval(nh_sim$attributable - s_sim$attributable)

  scenarios <- rbind(scenarios, data.frame(
    scenario         = s,
    episodes         = round(s_point$n_exposed),
    attributable     = round(s_point$attributable),
    prevented        = round(prevented),
    prevented_lower  = round(spread["lower"]),
    prevented_upper  = round(spread["upper"]),
    pct_attributable = round(100 * prevented / nh_point$attributable, 1),
    pct_all_asthma   = round(100 * prevented / nh_point$asthma_total, 1)))
}

cat("\nSCENARIOS AGAINST NATURAL HISTORY\n\n")
cat("natural history: episodes", round(nh_point$n_exposed),
    " attributable", round(nh_point$attributable), "\n\n")
print(scenarios, row.names = FALSE)


# One-way sensitivity on the base case odds ratio. Each input moves to
# the ends of its own interval while the others stay at base case.

odds_ratio <- effects$odds_ratio[1]
low_high   <- list(
  "odds ratio"        = c(effects$lower[1], effects$upper[1]),
  "RSV LRTI episodes" = unname(quantile(n_exposed_trials, c(0.025, 0.975))),
  "prevalence anchor" = c(prev_total_lower, prev_total_upper))

tornado <- data.frame()

for (input in names(low_high)) {

  or   <- c(odds_ratio, odds_ratio)
  nex  <- c(n_exposed_point, n_exposed_point)
  prev <- c(prev_total, prev_total)

  if (input == "odds ratio")        or   <- low_high[[input]]
  if (input == "RSV LRTI episodes") nex  <- low_high[[input]]
  if (input == "prevalence anchor") prev <- low_high[[input]]

  cal <- prevalence_unexposed(or, prev, nex / pop_total)
  out <- outcomes(nex, pop_total, cal$p0, cal$risk_ratio)

  tornado <- rbind(tornado, data.frame(
    input = input,
    low   = round(100 * out$attributable_fraction[1], 1),
    high  = round(100 * out$attributable_fraction[2], 1)))
}

tornado$width <- abs(tornado$high - tornado$low)
tornado <- tornado[order(-tornado$width), ]

cat("\nONE-WAY SENSITIVITY, attributable fraction in percent\n\n")
print(tornado, row.names = FALSE)


# Output

write.csv(results,   "results_effect_estimates.csv", row.names = FALSE)
write.csv(scenarios, "results_scenarios.csv",        row.names = FALSE)
write.csv(tornado,   "results_one_way.csv",          row.names = FALSE)
