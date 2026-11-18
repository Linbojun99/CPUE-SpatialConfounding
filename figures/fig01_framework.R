# ==============================================================================
# figures/fig01_framework.R
#
# Figure 1: diagram of the simulation framework. Five stages (conditioning,
# operating model, estimation, prediction, evaluation) are drawn as boxes
# joined by arrows, with a separate box for the Gaussian simulation. The
# diagram uses no data.
#
# Input:
#   none
# Output:
#   output/figures/fig01_framework.png
#   output/figures/fig01_framework.pdf
# ==============================================================================

source("00_settings.R")
source("figures/00_figure_style.R")
library(grid)


# 1. Size, coordinates and colours ---------------------------------------------

# Positions are in mm on a 170 x 92 mm page: x from the left edge and y from
# the top edge. The y coordinate of the page starts at the bottom, so every
# element is drawn at figure_height - y.
figure_height <- 92

fill_grey      <- "#E9EDEF"
fill_teal      <- "#DCEDE8"
fill_blue      <- "#D7E4F1"
fill_pale_blue <- "#ECF3FA"
fill_orange    <- "#F7E1CA"
colour_text    <- "#172C38"
colour_border  <- "#7A8991"
colour_arrow   <- "#39464D"

# Height of the coloured band that holds the title of a box.
band_height <- 6


# 2. Boxes ---------------------------------------------------------------------

# One row per box: top-left corner, width, height, colour of the title band,
# title, body text and font size of the body text. The body of the operating
# model box and the content of the prediction box are added in section 3.
boxes <- data.table(
  x         = c(6, 112, 6, 66, 119, 57),
  y         = c(2, 2, 24, 24, 24, 70),
  width     = c(46, 52, 50, 42, 45, 75),
  height    = c(14, 14, 30, 22, 28, 20),
  band_fill = c(fill_grey, fill_grey, fill_teal, fill_blue, fill_pale_blue,
                fill_orange),
  title     = c("1. Conditioning", "Gaussian simulation",
                "2. Operating model", "3. Estimation", "4. Prediction",
                "5. Evaluation"),
  body      = c(
    "Logbook data\nParameter estimates (Table 1)",
    "Four scenarios on a unit square\nSame stages 2-5 (Section 2.7)",
    "",
    paste("Six fitted models (Table 3)", "RSR transformation (Eq. 17)",
          "Convergence checks", sep = "\n"),
    "",
    paste("Coefficient error and interval coverage",
          "Index RMSE and change error",
          "Density error in cells with and without sets", sep = "\n")
  ),
  body_size = c(7.5, 7.2, 7.2, 7.5, 7.5, 7.2)
)

# Each box is a white rectangle with a border and, inside its upper edge, a
# coloured band without a border. The prediction box also holds two
# coloured panels for the two roles of a covariate.
rectangles <- rbind(
  boxes[, .(x, y, width, height, fill = "white", border = colour_border)],
  boxes[, .(x = x + 0.15, y = y + 0.15, width = width - 0.3,
            height = band_height, fill = band_fill, border = NA)],
  data.table(x = c(120.5, 142), y = 32.5, width = c(20, 20.5), height = 12.5,
             fill = fill_pale_blue, border = NA)
)


# 3. Text ----------------------------------------------------------------------

# One row per piece of text: position, font size in points, whether it is
# bold, and the horizontal justification (0 left, 0.5 centred). Titles are
# centred in the band and body text in the remaining area of the box.
scenario_groups <- c(
  A = "Confounding and collinearity",
  B = "Both covariates affect density",
  C = "No covariate effect",
  D = "Spatial scale of the covariate",
  E = "Trend and covariate role"
)
group_y <- 34.2 + 3.55 * seq_along(scenario_groups)

text_labels <- rbind(
  # Titles and body text of the boxes.
  boxes[, .(label = title, x = x + width / 2, y = y + 3.15, size = 8.5,
            bold = TRUE, hjust = 0.5)],
  boxes[body != "", .(label = body, x = x + width / 2,
                      y = y + 6.15 + (height - 6.15) / 2, size = body_size,
                      bold = FALSE, hjust = 0.5)],
  # Operating model: number of scenarios and what each group varies.
  data.table(label = "17 scenarios, 50 replicates each", x = 31, y = 34.2,
             size = 7.2, bold = FALSE, hjust = 0.5),
  data.table(label = names(scenario_groups), x = 12.2, y = group_y,
             size = 7.2, bold = FALSE, hjust = 0.5),
  data.table(label = scenario_groups, x = 15.6, y = group_y,
             size = 7.2, bold = FALSE, hjust = 0),
  # Prediction: treatment of density and catchability covariates.
  data.table(
    label = c("Density", "Retain values", "Catchability", "Fix at zero",
              "Integrate over the prediction domain"),
    x     = c(130.5, 130.5, 152.25, 152.25, 141.5),
    y     = c(36.2, 41.2, 36.2, 41.2, 48.9),
    size  = 7.5,
    bold  = c(TRUE, FALSE, TRUE, FALSE, FALSE),
    hjust = 0.5
  ),
  # Labels of the arrows.
  data.table(
    label = c("Catches", "Fits", "True values", "Coefficients",
              "Predictions"),
    x     = c(61, 113.5, 31, 89, 144),
    y     = c(32.6, 32.6, 76.5, 61, 67),
    size  = 7,
    bold  = FALSE,
    hjust = c(0.5, 0.5, 0, 0, 0)
  )
)


# 4. Arrows --------------------------------------------------------------------

# One row per point of an arrow; the arrow head is at the last point.
#   1 conditioning to operating model (parameter values)
#   2 operating model to estimation (catches)
#   3 estimation to prediction (fits)
#   4 operating model to evaluation (true values)
#   5 estimation to evaluation (coefficients)
#   6 prediction to evaluation (predictions)
arrow_points <- data.table(
  arrow = c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6),
  x     = c(29, 29, 56, 66, 108, 119, 29, 29, 57, 87, 87, 141.5, 141.5, 132),
  y     = c(16, 24, 35, 35, 35, 35, 54, 80, 80, 46, 70, 52, 80, 80)
)


# 5. Assemble the diagram ------------------------------------------------------

rectangle_grob <- rectGrob(
  x      = unit(rectangles$x, "native"),
  y      = unit(figure_height - rectangles$y, "native"),
  width  = unit(rectangles$width, "native"),
  height = unit(rectangles$height, "native"),
  just   = c("left", "top"),
  gp     = gpar(fill = rectangles$fill, col = rectangles$border, lwd = 0.65)
)

text_grob <- textGrob(
  label = text_labels$label,
  x     = unit(text_labels$x, "native"),
  y     = unit(figure_height - text_labels$y, "native"),
  hjust = text_labels$hjust,
  vjust = 0.5,
  gp    = gpar(fontfamily = font_family, fontsize = text_labels$size,
               fontface = ifelse(text_labels$bold, "bold", "plain"),
               col = colour_text, lineheight = 1.15)
)

arrow_grob <- polylineGrob(
  x     = unit(arrow_points$x, "native"),
  y     = unit(figure_height - arrow_points$y, "native"),
  id    = arrow_points$arrow,
  arrow = arrow(length = unit(1.25, "mm"), type = "closed"),
  gp    = gpar(col = colour_arrow, fill = colour_arrow, lwd = 0.75,
               linejoin = "round")
)

# Rectangles first, then the arrows, then the text.
figure <- gTree(
  children = gList(rectangle_grob, arrow_grob, text_grob),
  vp       = viewport(xscale = c(0, figure_width_double),
                      yscale = c(0, figure_height))
)


# 6. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figures, "fig01_framework.png"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = ragg::agg_png
)
ggsave(
  file.path(dir_figures, "fig01_framework.pdf"),
  plot   = figure,
  width  = figure_width_double,
  height = figure_height,
  units  = "mm",
  dpi    = figure_dpi,
  bg     = "white",
  device = cairo_pdf
)
