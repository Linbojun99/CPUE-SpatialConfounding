# ==============================================================================
# figures/00_figure_style.R
#
# Plotting settings shared by all figure scripts: figure sizes, fonts, line
# widths, colours, the colours, point shapes, order and labels of the
# estimators, labels of symbols, the themes and the settings of the maps.
# The file defines objects only and draws nothing.
#
# Every figure script starts with source("00_settings.R") followed by
# source("figures/00_figure_style.R") and is run with the working directory
# set to the folder that contains 00_settings.R, for example
#   Rscript figures/fig02_simulated_fields.R
# Each figure is written to output/figures/ as a PNG file (600 dpi) and a
# PDF file.
# ==============================================================================


# 1. Packages ------------------------------------------------------------------

library(ggplot2)
library(patchwork)


# 2. Figure size and resolution ------------------------------------------------

# Widths of a double-column and a single-column figure (mm) and the
# resolution of the PNG files (dots per inch).
figure_width_double <- 170
figure_width_single <- 85
figure_dpi          <- 600


# 3. Fonts ---------------------------------------------------------------------

# All text is set in Arial. Sizes are in points.
font_family           <- "Arial"
font_size_base        <- 9
font_size_axis_text   <- 7.5
font_size_axis_title  <- 8
font_size_panel_title <- 8.5
font_size_legend      <- 7.5

# geom_text() and annotate() take the text size in mm; a size in points is
# divided by this constant.
points_per_mm <- 2.845276


# 4. Lines ---------------------------------------------------------------------

# Panel borders, axis ticks and reference lines have the same width. Zero
# lines and lines of equality are solid grey; nominal levels, design values
# and true values are dashed black.
line_width_frame     <- 0.3
line_width_reference <- 0.3
colour_reference     <- "grey65"
linetype_dashed      <- "42"


# 5. Colours -------------------------------------------------------------------

# Simulated quantities that do not belong to an estimator, and a neutral grey.
colour_simulated <- "#4292C6"
colour_neutral   <- "#8C8C8C"

# Outline of boxes in boxplots.
colour_box_outline <- "#333333"

# End colours of the diverging scale of the maps (white at zero). The same
# two colours identify the two models fitted to the logbook data (Fig. S1).
colour_scale_low  <- "#2166AC"
colour_scale_high <- "#B2182B"


# 6. Estimators of the fishery simulation --------------------------------------

# Order in legends and panels. The codes are those of the column "model" in
# 00_settings.R; M7 is restricted spatial regression (RSR).
estimator_order <- c("M1", "M2", "M3", "M4", "M6", "M5", "M7")

# Labels as plotmath expressions, so that the subscripts of X1 and X2 are
# typeset. The empty subscript gives labels without a subscript the same
# height as the others.
estimator_labels <- c(
  M1 = "paste('GLMM (',X[1],', ',X[2],')')",
  M2 = "paste('ST (none)',''[phantom(1)])",
  M3 = "paste('ST (',X[1],')')",
  M4 = "paste('ST (',X[1],', ',X[2],')')",
  M6 = "paste('ST (',X[2],')')",
  M5 = "paste('ST spatial+',''[phantom(1)])",
  M7 = "paste('RSR',''[phantom(1)])"
)

# The same labels on two lines, for narrow columns.
estimator_labels_two_lines <- c(
  M1 = "atop('GLMM',paste('(',X[1],', ',X[2],')'))",
  M2 = "atop('ST','(none)')",
  M3 = "atop('ST',paste('(',X[1],')'))",
  M4 = "atop('ST',paste('(',X[1],', ',X[2],')'))",
  M6 = "atop('ST',paste('(',X[2],')'))",
  M5 = "atop('ST','spatial+')",
  M7 = "atop('RSR',phantom('(g)'))"
)

# Colours. In boxplots the boxes are filled with these colours and outlined
# with colour_box_outline.
estimator_colours <- c(
  M1 = "#8C8C8C",
  M2 = "#9ECAE1",
  M3 = "#4292C6",
  M4 = "#08306B",
  M6 = "#1A9641",
  M5 = "#F4A582",
  M7 = "#D6604D"
)

# Point shapes: open circle, filled circle, triangle, square, open diamond,
# inverted triangle and filled diamond. Shapes 21 to 25 are filled with the
# colour of the estimator.
estimator_shapes <- c(
  M1 = 1,
  M2 = 21,
  M3 = 24,
  M4 = 22,
  M6 = 5,
  M5 = 25,
  M7 = 23
)


# 7. Estimators of the Gaussian simulation -------------------------------------

# Codes of the column "model" in the output of 08a_gaussian_simulation.R, in
# the order of the panels, with the labels of the figure. Each estimator
# takes the colour and shape of its counterpart in the fishery simulation.
gaussian_estimator_order <- c("GLM", "SGLMM", "RSR", "spatial+")

gaussian_estimator_labels <- c(
  "GLM"      = "GLM",
  "SGLMM"    = "Spatial (TPRS)",
  "RSR"      = "RSR",
  "spatial+" = "spatial+"
)

gaussian_estimator_colours <- c(
  "GLM"      = "#8C8C8C",
  "SGLMM"    = "#08306B",
  "RSR"      = "#D6604D",
  "spatial+" = "#F4A582"
)

gaussian_estimator_shapes <- c(
  "GLM"      = 1,
  "SGLMM"    = 22,
  "RSR"      = 23,
  "spatial+" = 25
)


# 8. Symbols in labels ---------------------------------------------------------

# Greek letters are written as text inside plotmath expressions, so that they
# are set in Arial like the rest of the figure (the plotmath names rho, omega
# and beta would switch to a symbol font).
symbol_omega     <- "'\u03c9'"
symbol_rho_omega <- "'\u03c1'['\u03c9']"
symbol_rho_12    <- "'\u03c1'[12]"
symbol_beta_1    <- "'\u03b2'[1]"

# Negative numbers are printed with a minus sign (\u2212) instead of a
# hyphen. Axis labels get it from
#   scales::label_number(style_negative = "minus")
minus_sign <- "\u2212"


# 9. Axes ----------------------------------------------------------------------

# Year axes have ticks at the first, middle and last year. Narrow panels
# label the first and last year only.
year_breaks        <- c(2015, 2019, 2023)
year_labels_narrow <- c("2015", "", "2023")

# The 17 scenarios of Table 2, in the order of the table.
scenario_order <- main_scenarios


# 10. Themes -------------------------------------------------------------------

# Theme of all figures: black panel frame, no grid lines, panel titles in
# plain type on the left, legend below the panels without a title.
theme_figure <- theme_bw(
  base_size   = font_size_base,
  base_family = font_family
) +
  theme(
    panel.grid           = element_blank(),
    panel.border         = element_rect(colour = "black", fill = NA,
                                        linewidth = line_width_frame),
    axis.ticks           = element_line(colour = "black",
                                        linewidth = line_width_frame),
    axis.ticks.length    = unit(0.8, "mm"),
    axis.text            = element_text(size = font_size_axis_text,
                                        colour = "black"),
    axis.title           = element_text(size = font_size_axis_title),
    strip.background     = element_blank(),
    strip.text           = element_text(size = font_size_panel_title,
                                        face = "plain", hjust = 0),
    plot.title           = element_text(size = font_size_panel_title,
                                        face = "plain", hjust = 0),
    plot.tag             = element_text(size = font_size_panel_title,
                                        face = "plain"),
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.text          = element_text(size = font_size_legend),
    legend.key           = element_blank(),
    legend.key.size      = unit(3, "mm"),
    legend.spacing.x     = unit(1, "mm"),
    legend.key.spacing.x = unit(3, "mm"),
    legend.margin        = margin(0, 0, 0, 0),
    panel.spacing        = unit(2, "mm"),
    plot.margin          = margin(1, 1, 1, 1, "mm")
  )

# Theme of the maps: smaller axis text and a titled colour bar.
theme_figure_map <- theme_figure +
  theme(
    axis.text    = element_text(size = 7, colour = "black"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    strip.text   = element_text(size = font_size_panel_title, hjust = 0,
                                face = "plain", margin = margin(2, 1, 2, 1)),
    plot.margin  = margin(2, 2, 1, 2, "mm")
  )


# 11. Maps ---------------------------------------------------------------------

# Extent of the maps in degrees: the centres of the cells of the prediction
# domain span 59.5W to 13.5E and 32.5S to 28.5N, and the maps add 2 degrees
# on each side. Axis labels are placed at 50W, 25W and 0 and at 20S, 0 and
# 20N, inside the panels.
map_xlim       <- c(-61.5, 15.5)
map_ylim       <- c(-34.5, 30.5)
map_breaks_lon <- c(-50, -25, 0)
map_breaks_lat <- c(-20, 0, 20)

# Coastline: colour and line width.
colour_coastline     <- "grey50"
line_width_coastline <- 0.12

# Maps of standardized fields share a colour scale limited to -3 and 3;
# values beyond the limit take the end colours and the end labels carry
# the signs "at most" (\u2264) and "at least" (\u2265).
map_colour_limit  <- 3
map_colour_breaks <- -3:3
map_colour_labels <- c(
  paste0("\u2264 ", minus_sign, "3"),
  paste0(minus_sign, "2"),
  paste0(minus_sign, "1"),
  "0", "1", "2",
  "\u2265 3"
)

# Colour bar below the maps: title above the bar, bar 65 mm wide.
map_colour_bar <- guide_colourbar(
  title.position  = "top",
  title.hjust     = 0.5,
  barwidth        = unit(65, "mm"),
  barheight       = unit(2.5, "mm"),
  ticks.colour    = "#333333",
  frame.colour    = "#777777",
  frame.linewidth = 0.2
)
