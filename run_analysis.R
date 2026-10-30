# Edit the assumptions, then run the whole file.

source("rsv_asthma_model.R")


# Assumptions

set.seed(20261017)
trials <- 10000

input_cache <- "inputs.rds"

pop_total <- 3664292          # NCHS 2021 births

prev_total       <- 0.088     # NHIS 2018-2021, asthma, ages 5-11
prev_total_lower <- 0.081
prev_total_upper <- 0.096

episode_private  <- 69068     # Gantenberg, visit to episode correction
visit_private    <- 107338
episode_medicaid <- 170744
visit_medicaid   <- 256495
share_private    <- 0.491
share_medicaid   <- 0.374

# Narrow case definition, Gantenberg, specific ICD-10 codes
episode_private_narrow  <- 25409
visit_private_narrow    <- 43301
episode_medicaid_narrow <- 67357
visit_medicaid_narrow   <- 110167

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
encounters <- inputs$scenarios[["Natural History"]]

for (s in names(inputs$scenarios)) {
  if (length(inputs$scenarios[[s]]$trials) != trials)
    stop("Expected ", trials, " trials in ", input_cache)
}

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

# Normal approximation to the reported prevalence interval.
prev_total_trials <- rnorm(trials, prev_total,
                           (prev_total_upper - prev_total_lower) / (2 * 1.96))

# One set of normal draws, reused for every odds ratio, so rows differ
# only by the effect estimate.
z <- rnorm(trials)

# Narrow case definition. These draws come last so the earlier draws
# are unchanged.
correction_point_narrow <- w_private  * episode_private_narrow  / visit_private_narrow +
                           w_medicaid * episode_medicaid_narrow / visit_medicaid_narrow

correction_trials_narrow <-
  w_private  * rbeta(trials, episode_private_narrow + 1,
                     visit_private_narrow - episode_private_narrow + 1) +
  w_medicaid * rbeta(trials, episode_medicaid_narrow + 1,
                     visit_medicaid_narrow - episode_medicaid_narrow + 1)


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


# Table 1. Base case odds ratio. Point estimate from the base case,
# interval from the trials.

table1_draws <- cbind(
  encounters$by_setting[["Outpatient"]]$trials,
  encounters$by_setting[["Emergency Department"]]$trials,
  encounters$by_setting[["Hospitalizations"]]$trials,
  encounters$trials,
  nh_sim$n_exposed,
  100 * nh_sim$p_exposed,
  100 * nh_sim$prev_unexposed,
  100 * nh_sim$prev_exposed,
  nh_sim$risk_ratio,
  nh_sim$asthma_total,
  nh_sim$asthma_if_none_exposed,
  nh_sim$attributable,
  100000 * nh_sim$attributable / pop_total,
  100 * nh_sim$attributable_fraction)

table1 <- data.frame(
  quantity = c("encounters, outpatient", "encounters, emergency department",
               "encounters, hospitalization", "encounters, all settings",
               "unique illness episodes", "incidence, % of birth cohort",
               "asthma prevalence, no RSV LRTI (P0), %",
               "asthma prevalence, RSV LRTI (P1), %",
               "relative risk",
               "all-cause cases", "all-cause cases, no RSV LRTI",
               "RSV-attributable cases",
               "RSV-attributable cases per 100,000",
               "attributable fraction, %"),
  point = c(unname(encounters$by_setting[["Outpatient"]]$point),
            unname(encounters$by_setting[["Emergency Department"]]$point),
            unname(encounters$by_setting[["Hospitalizations"]]$point),
            unname(encounters$point),
            nh_point$n_exposed,
            100 * nh_point$p_exposed,
            100 * nh_point$prev_unexposed,
            100 * nh_point$prev_exposed,
            nh_point$risk_ratio,
            nh_point$asthma_total,
            nh_point$asthma_if_none_exposed,
            nh_point$attributable,
            100000 * nh_point$attributable / pop_total,
            100 * nh_point$attributable_fraction),
  lower = 0,
  upper = 0)

for (i in seq_len(nrow(table1))) {
  table1$lower[i] <- quantile(table1_draws[, i], 0.025)
  table1$upper[i] <- quantile(table1_draws[, i], 0.975)
}

cat("\nTABLE 1\n\n")
print(data.frame(quantity = table1$quantity,
                 point = round(table1$point, 2),
                 lower = round(table1$lower, 2),
                 upper = round(table1$upper, 2)), row.names = FALSE)


# Table 2. Nirsevimab at 50% coverage against natural history.
# Calibrated once on natural history and held fixed.

nir_point <- outcomes(unname(inputs$scenarios[["Nirsevimab"]]$point) * correction_point,
                      pop_total, cal_point$p0, cal_point$risk_ratio)
nir_sim   <- outcomes(unname(inputs$scenarios[["Nirsevimab"]]$trials) * correction_trials,
                      pop_total, cal_sim$p0, cal_sim$risk_ratio)

table2_draws <- cbind(
  nh_sim$n_exposed,       nir_sim$n_exposed,
  100 * (1 - nir_sim$n_exposed / nh_sim$n_exposed),
  nh_sim$asthma_total,    nir_sim$asthma_total,
  100 * (1 - nir_sim$asthma_total / nh_sim$asthma_total),
  nh_sim$attributable,    nir_sim$attributable,
  100 * (1 - nir_sim$attributable / nh_sim$attributable),
  nh_sim$attributable - nir_sim$attributable)

table2 <- data.frame(
  quantity = c("episodes, natural history", "episodes, nirsevimab",
               "episodes, decrease %",
               "all-cause cases, natural history", "all-cause cases, nirsevimab",
               "all-cause cases, decrease %",
               "attributable cases, natural history",
               "attributable cases, nirsevimab",
               "attributable cases, decrease %",
               "cases prevented"),
  point = c(nh_point$n_exposed, nir_point$n_exposed,
            100 * (1 - nir_point$n_exposed / nh_point$n_exposed),
            nh_point$asthma_total, nir_point$asthma_total,
            100 * (1 - nir_point$asthma_total / nh_point$asthma_total),
            nh_point$attributable, nir_point$attributable,
            100 * (1 - nir_point$attributable / nh_point$attributable),
            nh_point$attributable - nir_point$attributable),
  lower = 0,
  upper = 0)

for (i in seq_len(nrow(table2))) {
  table2$lower[i] <- quantile(table2_draws[, i], 0.025)
  table2$upper[i] <- quantile(table2_draws[, i], 0.975)
}

cat("\nTABLE 2\n\n")
print(data.frame(quantity = table2$quantity,
                 point = round(table2$point, 1),
                 lower = round(table2$lower, 1),
                 upper = round(table2$upper, 1)), row.names = FALSE)


# Other coverage and effectiveness assumptions. Each workbook is
# calibrated on its own natural history trials and compared with its
# own nirsevimab trials.

arms <- c("90% coverage"                 = "inputs_90pct.rds",
          "post-licensure effectiveness" = "inputs_postlicensure.rds")

other_arms <- data.frame()

for (arm in names(arms)) {

  arm_inputs <- readRDS(arms[[arm]])$scenarios

  arm_nh_point  <- unname(arm_inputs[["Natural History"]]$point)  * correction_point
  arm_nh_trials <- arm_inputs[["Natural History"]]$trials         * correction_trials

  arm_cal_point <- prevalence_unexposed(effects$odds_ratio[1], prev_total,
                                        arm_nh_point / pop_total)
  arm_cal_sim   <- prevalence_unexposed(odds_draws, prev_total_trials,
                                        arm_nh_trials / pop_total)

  arm_nh_out  <- outcomes(arm_nh_point, pop_total,
                          arm_cal_point$p0, arm_cal_point$risk_ratio)
  arm_nh_sim  <- outcomes(arm_nh_trials, pop_total,
                          arm_cal_sim$p0, arm_cal_sim$risk_ratio)

  arm_nir_out <- outcomes(unname(arm_inputs[["Nirsevimab"]]$point) * correction_point,
                          pop_total, arm_cal_point$p0, arm_cal_point$risk_ratio)
  arm_nir_sim <- outcomes(arm_inputs[["Nirsevimab"]]$trials * correction_trials,
                          pop_total, arm_cal_sim$p0, arm_cal_sim$risk_ratio)

  prevented     <- arm_nh_out$attributable - arm_nir_out$attributable
  prevented_sim <- arm_nh_sim$attributable - arm_nir_sim$attributable
  pct_attrib    <- 100 * prevented_sim / arm_nh_sim$attributable
  pct_all       <- 100 * prevented_sim / arm_nh_sim$asthma_total

  other_arms <- rbind(other_arms, data.frame(
    arm              = arm,
    prevented        = round(prevented),
    prevented_lower  = round(quantile(prevented_sim, 0.025)),
    prevented_upper  = round(quantile(prevented_sim, 0.975)),
    pct_attributable = round(100 * prevented / arm_nh_out$attributable, 1),
    pct_attrib_lower = round(quantile(pct_attrib, 0.025), 1),
    pct_attrib_upper = round(quantile(pct_attrib, 0.975), 1),
    pct_all_asthma   = round(100 * prevented / arm_nh_out$asthma_total, 1),
    pct_all_lower    = round(quantile(pct_all, 0.025), 1),
    pct_all_upper    = round(quantile(pct_all, 0.975), 1)))
}

cat("\nOTHER COVERAGE AND EFFECTIVENESS ASSUMPTIONS\n\n")
print(other_arms, row.names = FALSE)


# Narrow case definition. Only the encounter-to-episode correction
# changes.

narrow_episodes  <- unname(encounters$point) * correction_point_narrow
cal_narrow_point <- prevalence_unexposed(effects$odds_ratio[1], prev_total,
                                         narrow_episodes / pop_total)
narrow_base <- outcomes(narrow_episodes, pop_total,
                        cal_narrow_point$p0, cal_narrow_point$risk_ratio)

narrow_trials   <- encounters$trials * correction_trials_narrow
cal_narrow_sim  <- prevalence_unexposed(odds_draws, prev_total_trials,
                                        narrow_trials / pop_total)
narrow_sim <- outcomes(narrow_trials, pop_total,
                       cal_narrow_sim$p0, cal_narrow_sim$risk_ratio)

narrow_fraction <- interval(100 * narrow_sim$attributable_fraction)

# Attributable cases under the narrow definition. The Results report
# this count beside the attributable fraction.
narrow_cases <- interval(narrow_sim$attributable)

cat("\nNARROW CASE DEFINITION\n\n")
cat("correction factor", round(100 * correction_point_narrow, 1), "%\n")
cat("attributable fraction", round(100 * narrow_base$attributable_fraction, 1),
    "% (", round(narrow_fraction["lower"], 1), "to",
    round(narrow_fraction["upper"], 1), ")\n")
cat("attributable cases", round(narrow_base$attributable),
    "(", round(narrow_cases["lower"]), "to",
    round(narrow_cases["upper"]), ")\n")

narrow <- data.frame(
  correction_factor = round(100 * correction_point_narrow, 1),
  fraction          = round(100 * narrow_base$attributable_fraction, 2),
  fraction_lower    = round(narrow_fraction["lower"], 1),
  fraction_upper    = round(narrow_fraction["upper"], 1),
  cases             = round(narrow_base$attributable),
  cases_lower       = round(narrow_cases["lower"]),
  cases_upper       = round(narrow_cases["upper"]))


# One-way sensitivity on the base case odds ratio. Each input moves to
# the ends of its own interval while the others stay at base case.

odds_ratio <- effects$odds_ratio[1]
low_high   <- list(
  "odds ratio"        = c(effects$lower[1], effects$upper[1]),
  "RSV LRTI encounters" =
    correction_point * unname(quantile(encounters$trials, c(0.025, 0.975))),
  "prevalence anchor" = c(prev_total_lower, prev_total_upper),
  "encounter-to-episode correction" =
    unname(encounters$point) * unname(quantile(correction_trials, c(0.025, 0.975))))

tornado <- data.frame()

for (input in names(low_high)) {

  or   <- c(odds_ratio, odds_ratio)
  nex  <- c(n_exposed_point, n_exposed_point)
  prev <- c(prev_total, prev_total)

  if (input == "odds ratio")        or   <- low_high[[input]]
  if (input == "RSV LRTI encounters") nex  <- low_high[[input]]
  if (input == "prevalence anchor") prev <- low_high[[input]]
  if (input == "encounter-to-episode correction") nex <- low_high[[input]]

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
write.csv(table1,    "results_table1.csv",           row.names = FALSE)
write.csv(table2,    "results_table2.csv",           row.names = FALSE)
write.csv(other_arms, "results_other_arms.csv",      row.names = FALSE)
write.csv(narrow,     "results_narrow.csv",          row.names = FALSE)
