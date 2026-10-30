# Figure 1: all-cause asthma at age 6 (panel A) and RSV LRTI-attributable
# asthma split into cases prevented by nirsevimab at 50% coverage and
# residual attributable cases not prevented (panel B).
#
# Every number on the figure is read from the results files that
# run_analysis.R writes, so the figure changes when the analysis changes.
# Nothing is typed in by hand. Run run_analysis.R first, then this file,
# from the same working folder. Base R only, no extra packages.

# ---- Settings ----------------------------------------------------------

table1_file <- "results_table1.csv"
table2_file <- "results_table2.csv"
out_pdf     <- "Figure1.pdf"
out_png     <- "Figure1.png"

col_not_attrib <- "grey95"    # asthma not attributable to RSV LRTI
# Light fills so black text reads on every bar.
col_prevented  <- "#F4B183"   # attributable cases prevented by nirsevimab
col_residual   <- "#9DC3E6"   # attributable cases not prevented

top_a <- 350000   # top of the panel A axis
top_b <- 50000    # top of the panel B axis

# Panel B sits higher than panel A so a reader can see that the two
# panels do not share a vertical axis. Both settings are shares of the
# panel A height. b_base is where the panel B zero line sits, and
# b_height is how tall the panel B axis is, so the top of bar B does not
# line up with the top of bar A.
b_base   <- 0.25
b_height <- 0.60

fig_width  <- 10.9 # inches
fig_height <- 4.5  # inches

# Arial matches the manuscript. The cairo devices used below embed the
# font in the PDF.
font_family <- "Arial"

# ---- Read the results --------------------------------------------------

t1 <- read.csv(table1_file)
t2 <- read.csv(table2_file)

# Pull one row (point, lower, upper) by its quantity name.
get_row <- function(tab, name) {
  row <- tab[tab$quantity == name, c("point", "lower", "upper")]
  if (nrow(row) != 1) stop("Not found in results file: ", name)
  unlist(row)
}

all_cause   <- get_row(t1, "all-cause cases")
attributable <- get_row(t1, "RSV-attributable cases")
paf         <- get_row(t1, "attributable fraction, %")
prevented   <- get_row(t2, "cases prevented")
not_prev    <- get_row(t2, "attributable cases, nirsevimab")
pct_prev    <- get_row(t2, "attributable cases, decrease %")

# The share not prevented is 100 minus the share prevented. Its interval
# limits swap, because the lower limit of one is the upper limit of the other.
pct_not_prev <- c(100 - pct_prev[1], 100 - pct_prev[3], 100 - pct_prev[2])

# ---- Formatting --------------------------------------------------------

# Counts to whole numbers with commas. Percentages to one decimal with a
# middle dot as the decimal point, as the Lancet journals print them.
# Rounding happens only here, from the unrounded values in the files.
fmt_n <- function(x) format(round(x), big.mark = ",", trim = TRUE)
fmt_p <- function(x) gsub(".", "\u00b7", sprintf("%.1f", x), fixed = TRUE)
ui_n  <- function(v) paste0(fmt_n(v[1]), " (", fmt_n(v[2]), "\u2013", fmt_n(v[3]), ")")
ui_p  <- function(v) paste0(fmt_p(v[1]), "% (", fmt_p(v[2]), "\u2013", fmt_p(v[3]), ")")

# ---- Drawing -----------------------------------------------------------

# The figure is drawn on one blank canvas with a vertical scale from 0
# to 1. Panel A maps 0-350,000 cases onto the full height. Panel B maps
# 0-50,000 cases onto a shorter space that starts at b_base, so its zero
# line sits above the zero line of panel A.
ya <- function(v) v / top_a
yb <- function(v) b_base + b_height * v / top_b

# A bracket from y0 to y1 at x, with a point at its middle facing right.
bracket <- function(x, y0, y1, w = 0.6) {
  segments(x, y0, x, y1, lwd = 1.5)
  segments(x - w, y0, x, y0, lwd = 1.5)
  segments(x - w, y1, x, y1, lwd = 1.5)
  segments(x, (y0 + y1) / 2, x + w, (y0 + y1) / 2, lwd = 1.5)
}

draw_figure <- function() {
  par(mar = c(0, 0, 0, 0), family = font_family)
  plot.new()
  # The canvas ends at x = 91, just past the right-hand labels. Changing
  # this needs a matching change to fig_width so the drawing keeps its shape.
  plot.window(xlim = c(0, 91), ylim = c(-0.16, 1.1))
  cx <- 0.72   # text size

  # Panel A: one stacked bar for all-cause asthma.
  xa <- c(12, 26)
  rect(xa[1], 0, xa[2], ya(not_prev[1]), col = col_residual)
  rect(xa[1], ya(not_prev[1]), xa[2], ya(attributable[1]), col = col_prevented)
  rect(xa[1], ya(attributable[1]), xa[2], ya(all_cause[1]), col = col_not_attrib)
  ticks_a <- seq(0, top_a, 50000)
  text(9, ya(ticks_a), fmt_n(ticks_a), adj = 1, cex = cx)
  text(3, ya(top_a / 2), "Asthma cases", srt = 90, cex = cx)
  text(mean(xa), ya((attributable[1] + all_cause[1]) / 2),
       "Asthma not\nattributable to\nRSV LRTI", cex = cx)
  text(mean(xa), -0.06, "All-cause cases", cex = cx)
  text(mean(xa), -0.11, ui_n(all_cause), cex = cx)
  # Both panel titles sit at the same height, each starting at the left
  # edge of its bar.
  text(xa[1], 1.06, "A. All-cause asthma at age 6", adj = 0, font = 2, cex = 0.85)

  # Panel B: the attributable cases on their own scale.
  xb <- c(38, 52)
  rect(xb[1], yb(0), xb[2], yb(not_prev[1]), col = col_residual)
  rect(xb[1], yb(not_prev[1]), xb[2], yb(attributable[1]), col = col_prevented)
  ticks_b <- seq(0, top_b, 10000)
  text(57, yb(ticks_b), fmt_n(ticks_b), adj = 1, cex = cx)
  # Axis title to the right of the panel B numbers, in the space before
  # the brackets.
  text(58.3, yb(top_b / 2), "Asthma cases", srt = 270, cex = cx)
  text(mean(xb), yb(not_prev[1] / 2), "Residual\nattributable\nasthma",
       cex = cx)
  text(mean(xb), yb((not_prev[1] + attributable[1]) / 2),
       "RSV-attributable\nasthma cases\nprevented", cex = cx)
  text(mean(xb), yb(0) - 0.06, "RSV LRTI\nattributable cases", cex = cx)
  text(mean(xb), yb(0) - 0.13, ui_n(attributable), cex = cx)
  text(xb[1], 1.06, "B. RSV LRTI-attributable asthma", adj = 0, font = 2, cex = 0.85)
  # The coverage note lines up with the start of the title text after "B. ".
  b_indent <- strwidth("B. ", font = 2, cex = 0.85)
  text(xb[1] + b_indent, 0.99, "(Nirsevimab at 50% coverage)", adj = 0, cex = cx)

  # Dashed lines that link the attributable part of bar A to bar B.
  segments(xa[2], ya(attributable[1]), xb[1], yb(attributable[1]), lty = 2)
  segments(xa[2], 0, xb[1], yb(0), lty = 2)

  # Brackets and labels to the right of panel B. Each label starts just
  # past the point of its bracket, and the lines of a label are kept
  # together so they sit beside the bar segment they describe.
  bracket(61.5, yb(not_prev[1]) + 0.01, yb(attributable[1]))
  text(62.8, yb((not_prev[1] + attributable[1]) / 2), adj = 0,
       paste0("Asthma cases prevented\n", ui_n(prevented),
              "\nAttributable burden prevented\n", ui_p(pct_prev)),
       cex = cx)

  bracket(61.5, yb(0), yb(not_prev[1]) - 0.01)
  text(62.8, yb(not_prev[1] / 2), adj = 0,
       paste0("Residual attributable cases\n", ui_n(not_prev),
              "\nAttributable burden remaining\n", ui_p(pct_not_prev)),
       cex = cx)

  bracket(81, yb(0), yb(attributable[1]))
  text(82.3, yb(attributable[1] / 2), adj = 0,
       paste0("Population\nattributable fraction\n", ui_p(paf)), cex = cx)
}

# ---- Save --------------------------------------------------------------

# cairo_pdf embeds the font and prints the en dash and the middle dot.
cairo_pdf(out_pdf, width = fig_width, height = fig_height, family = font_family)
draw_figure()
dev.off()

png(out_png, width = fig_width, height = fig_height, units = "in", res = 300,
    type = "cairo", family = font_family)
draw_figure()
dev.off()
