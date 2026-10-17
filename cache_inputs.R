# Reads one incidence workbook and saves the first-year encounter
# counts the asthma model needs. Run once per workbook.
#
# Trials are matched across scenarios inside a workbook. They are not
# matched across workbooks. Never subtract results from two workbooks.

library(readxl)

xlsx     <- "healthoutcomesnewsamples.xlsx"
out_file <- "inputs.rds"

scenarios <- c("Natural History", "Nirsevimab", "RSVpreF")
settings  <- c("Outpatient", "Emergency Department", "Hospitalizations")
months    <- as.character(0:11)

raw <- as.data.frame(read_excel(xlsx))

first_year <- function(label) rowSums(raw[raw$Intervention == label, months])

cached <- list()

for (scenario in scenarios) {

  trials     <- 0
  point      <- 0
  by_setting <- list()

  for (setting in settings) {
    label    <- paste(scenario, setting)
    s_trials <- first_year(label)
    s_point  <- first_year(paste(label, "PE"))
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
