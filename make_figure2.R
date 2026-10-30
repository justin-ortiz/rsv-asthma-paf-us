# Figure 2: one-way sensitivity analysis of the population attributable
# fraction. Each bar spans the attributable fraction when one input is at
# the low and high ends of its interval, with the others at base case.
#
# The bar ends come from results_one_way.csv and the base case line from
# results_table1.csv, both written by run_analysis.R. Nothing is typed in
# by hand. Run run_analysis.R first, then this file, from the same
# working folder. Base R only, no extra packages.

# ---- Settings ----------------------------------------------------------

one_way_file <- "results_one_way.csv"
table1_file  <- "results_table1.csv"
out_pdf      <- "Figure2.pdf"
out_png      <- "Figure2.png"

# Bar fill. Mid grey keeps the narrowest bar visible beside the base
# case line.
col_bar <- "#A6A6A6"

fig_width  <- 9     # inches
fig_height <- 3.6   # inches

# Labels printed on the figure for each input name in results_one_way.csv.
labels <- c("odds ratio"                      = "Adjusted odds ratio",
            "RSV LRTI encounters"             = "RSV LRTI incidence",
            "prevalence anchor"               = "Asthma prevalence anchor",
            "encounter-to-episode correction" = "Encounter-to-episode correction")

# ---- Read the results --------------------------------------------------

one_way <- read.csv(one_way_file)
t1      <- read.csv(table1_file)

base <- t1$point[t1$quantity == "attributable fraction, %"]
if (length(base) != 1) stop("Not found in results file: attributable fraction, %")

# The low end of an input does not always give the lower fraction. A higher
# prevalence anchor gives a lower fraction, for example. So each bar runs
# from the smaller to the larger of the two values.
one_way$left  <- pmin(one_way$low, one_way$high)
one_way$right <- pmax(one_way$low, one_way$high)

# run_analysis.R already sorts the rows from widest to narrowest. The
# widest bar goes at the top.
n <- nrow(one_way)
one_way$y <- n:1

# ---- Formatting --------------------------------------------------------

# One decimal with a middle dot as the decimal point, as the Lancet
# journals print it.
fmt_p <- function(x) gsub(".", "\u00b7", sprintf("%.1f", x), fixed = TRUE)

# ---- Drawing -----------------------------------------------------------

draw_figure <- function() {
  par(mar = c(4, 15, 2.5, 8), family = "sans", cex = 0.85)
  plot(NA, xlim = c(0, 25), ylim = c(0.5, n + 0.5),
       xaxs = "i", yaxt = "n", xaxt = "n", bty = "n",
       xlab = "Population attributable fraction, %", ylab = "")

  # Light vertical grid lines at each tick.
  ticks <- seq(0, 25, 5)
  abline(v = ticks, col = "grey88")
  axis(1, at = ticks, labels = ticks, lwd = 0, lwd.ticks = 0)

  rect(one_way$left, one_way$y - 0.3, one_way$right, one_way$y + 0.3,
       col = col_bar, border = NA)

  axis(2, at = one_way$y, labels = labels[one_way$input],
       las = 1, lwd = 0, lwd.ticks = 0)

  # Base case line, drawn last so it sits on top of the bars.
  segments(base, 0.5, base, n + 0.5, lty = 2, lwd = 1.5)
  text(base, n + 0.5, paste0("Base case ", fmt_p(base), "%"),
       pos = 3, xpd = NA)

  # Right-hand column with each bar's end values, so the two narrow bars
  # can be read without changing the scale.
  x_col <- 26.5
  text(x_col, n + 0.5, "Range %", adj = c(0, 0), xpd = NA)
  text(x_col, one_way$y,
       paste0(fmt_p(one_way$left), "\u2013", fmt_p(one_way$right)),
       adj = 0, xpd = NA)
}

# ---- Save --------------------------------------------------------------

# The WinAnsi encoding lets the PDF print the middle dot.
pdf(out_pdf, width = fig_width, height = fig_height, encoding = "WinAnsi.enc")
draw_figure()
dev.off()

png(out_png, width = fig_width, height = fig_height, units = "in", res = 300)
draw_figure()
dev.off()
