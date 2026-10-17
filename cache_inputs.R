# Reads one incidence workbook and saves the first-year encounter
# counts the asthma model needs. Run once per workbook.
#
# Trials are matched across scenarios inside a workbook, not across
# workbooks.
#
# Layout of the workbook, first sheet: six header rows, then the base
# case row, then one row per trial. Each block of 36 columns holds
# months 0 to 35 for one scenario and care setting. The numbers below
# are the first column of each LRTI block.

library(readxl)

xlsx     <- "RSV_incidence_50pct_coverage_base_waning.xlsx"
out_file <- "inputs.rds"

starts <- list(
  "Natural History" = c("Outpatient" = 76, "Emergency Department" = 112,
                        "Hospitalizations" = 148),
  "RSVpreF"         = c("Outpatient" = 256, "Emergency Department" = 292,
                        "Hospitalizations" = 328),
  "Nirsevimab"      = c("Outpatient" = 436, "Emergency Department" = 472,
                        "Hospitalizations" = 508))

raw <- suppressMessages(suppressWarnings(as.data.frame(
  read_excel(xlsx, sheet = "Monte Carlo Simulation Results",
             col_names = FALSE, skip = 6, col_types = "numeric"))))

if (nrow(raw) != 10001) stop("Expected the base case row and 10,000 trials.")

# Months 0 to 11 are the first year of life.
first_year <- function(start) rowSums(raw[, start:(start + 11)])

cached <- list()

for (scenario in names(starts)) {

  trials     <- 0
  point      <- 0
  by_setting <- list()

  for (setting in names(starts[[scenario]])) {
    s        <- first_year(starts[[scenario]][[setting]])
    s_point  <- s[1]
    s_trials <- s[-1]
    by_setting[[setting]] <- list(trials = s_trials, point = s_point)
    trials <- trials + s_trials
    point  <- point  + s_point
  }

  cached[[scenario]] <- list(trials = trials, point = point,
                             by_setting = by_setting)
}

saveRDS(list(source = basename(xlsx), scenarios = cached), out_file)

cat("cached", out_file, "with",
    length(cached[["Natural History"]]$trials), "trials\n")
cat("natural history, first-year encounters, base case\n")
for (setting in names(starts[["Natural History"]])) {
  cat(" ", setting,
      round(cached[["Natural History"]]$by_setting[[setting]]$point), "\n")
}
cat("  total", round(cached[["Natural History"]]$point), "\n")
