# One figure per case showing only how a find's date gets recorded, with the
# measured value left out. Every figure has the same layout: a top strip with the
# dating scheme the case uses (none for Case 2, which has no scheme), then a
# handful of simulated finds, one per row, sorted by true date. The y axis is
# only that ordering, not the measured value.
#
# Parameters are picked for legibility, not taken from the design grid. Cases 3
# and 4 share one uneven periodisation (alpha_conc = 1, seed 35) so the two
# figures differ only in how a window is read off it.

library(here)
suppressMessages(library(rcarbon))
library(ggplot2)
library(patchwork)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))
source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))
source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))
source(here("Simulations", "Sim_Case4_MergedPhases", "scripts", "simulate.R"))

OUT_DIR  <- here("Simulations", "shared", "figures")
ACCENT   <- "#780000"
N_SHOW   <- 16
TRUE_LAB <- "true date (simulated)"
ROW_LAB  <- "finds, sorted by true date"

base_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        legend.position = "bottom", legend.title = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

strip_theme <- base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

save_fig <- function(strip, rows, title, subtitle, file, xlim, strip_h = 1) {
  rows <- rows + coord_cartesian(xlim = xlim)
  fig  <- if (is.null(strip)) rows else
    (strip + coord_cartesian(xlim = xlim)) / rows +
      plot_layout(heights = c(strip_h, 3.2), guides = "collect")
  fig  <- fig + plot_annotation(title = title, subtitle = subtitle,
                                theme = base_theme)
  ggsave(file.path(OUT_DIR, file), fig & theme(legend.position = "bottom"),
         width = 9, height = if (is.null(strip)) 5.4 else 6.1 + 0.7 * strip_h, dpi = 300,
         bg = "white")
  cat("wrote", file, "\n")
}

sort_rows <- function(d) {
  d <- d[order(d$True_date), ]
  d$row <- seq_len(nrow(d))
  d
}

bce_axis <- scale_x_continuous(
  labels = function(x) ifelse(x < 0, paste0(abs(x), " BCE"), paste(x, "CE")))

# Case 1: calibrated distributions under the IntCal20 curve

window <- c(-800, -400)
c14 <- sort_rows(simulate_c14(N_SHOW, 8, 0.02, 2, window, lab_error = 30,
                              seed = 11))
cal <- calmatrix_to_calendar(calibrate(c14$CRA, errors = c14$Error,
                                       calCurves = "intcal20", calMatrix = TRUE,
                                       verbose = FALSE))
keep <- cal$years >= -1000 & cal$years <= -200
dens <- do.call(rbind, lapply(seq_len(nrow(c14)), function(i) {
  p <- cal$prob[keep, i]
  data.frame(row = i, year = cal$years[keep], ymin = i,
             ymax = i + 0.85 * p / max(p))
}))

yrs   <- seq(-1000, -200, by = 5)
cc    <- uncalibrate(1950 - yrs, verbose = FALSE)
curve <- data.frame(year = yrs, cra = cc$ccCRA, err = cc$ccError)

strip1 <- ggplot(curve, aes(year, cra)) +
  geom_rect(aes(xmin = window[1], xmax = window[2], ymin = -Inf, ymax = Inf,
                fill = "study window"), data = data.frame(x = 1), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = cra - 2 * err, ymax = cra + 2 * err), fill = "grey45",
              alpha = 0.6) +
  geom_line(linewidth = 0.4) +
  # Labelled on the curve, so its swatch does not sit next to the calibrated
  # distributions in the legend
  annotate("text", x = -250, y = max(curve$cra), label = "IntCal20, mean ± 2σ",
           hjust = 1, vjust = 1, size = 3.2, colour = "grey20") +
  scale_fill_manual(values = c("study window" = "grey92")) +
  labs(x = NULL, y = expression({}^14*C~age~(BP))) +
  strip_theme +
  theme(axis.text.y = element_text(), axis.ticks.y = element_line(),
        axis.line.y = element_line())

rows1 <- ggplot(dens) +
  geom_ribbon(aes(year, ymin = ymin, ymax = ymax, group = row,
                  fill = "calibrated date"),
              colour = "grey40", linewidth = 0.2) +
  geom_point(data = c14, aes(True_date, row, colour = TRUE_LAB), size = 2) +
  scale_fill_manual(values = c("calibrated date" = "grey70")) +
  scale_colour_manual(values = ACCENT) +
  bce_axis +
  labs(x = "calendar year", y = ROW_LAB) +
  base_theme

save_fig(strip1, rows1, "Case 1: radiocarbon",
         "Hallstatt plateau, lab error 30 years",
         "case1_dating.png", c(-1000, -200))

# Case 2: independent windows, no shared scheme

typo <- sort_rows(simulate_typo(N_SHOW, 8, 0.02, 2, prop_coarse_samples = 0.5, seed = 5))

typo$window <- ifelse(typo$Coarse, "coarsely dated window", "well-dated window")

rows2 <- ggplot(typo) +
  geom_segment(aes(x = Start_date, xend = End_date, y = row, yend = row,
                   colour = window), linewidth = 2.2) +
  geom_point(aes(True_date, row, colour = TRUE_LAB), size = 2) +
  scale_colour_manual(values = c("well-dated window" = "grey35",
                                 "coarsely dated window" = "grey75",
                                 setNames(ACCENT, TRUE_LAB)),
                      breaks = c("well-dated window", "coarsely dated window", TRUE_LAB)) +
  labs(x = "calendar year", y = ROW_LAB) +
  base_theme

save_fig(NULL, rows2, "Case 2: typochronology",
         "Each find has its own window: 6% of the period if well dated, 20-30% if coarsely dated",
         "case2_dating.png", c(-300, 1300))

K <- 6; ALPHA <- 1; PHASE_SEED <- 35

# Case 3: overlapping phases

ov    <- simulate_overlap(N_SHOW + 4, 8, 0.02, 2, K, ALPHA, overlap = 0.5,
                          assign_p = 0.5, seed = PHASE_SEED)
b3    <- attr(ov, "bounds")
reach <- 0.5 * diff(b3) / 2
ph3   <- data.frame(phase = seq_len(K), lo = b3[-(K + 1)] - reach,
                    hi = b3[-1] + reach, mid = (b3[-(K + 1)] + b3[-1]) / 2)
ov    <- sort_rows(ov)
in_overlap <- sapply(ov$True_date, function(t)
  any(t >= ph3$lo[-1] & t <= ph3$hi[-K]))
ov$where <- ifelse(in_overlap, "true date where two phases overlap",
                   "true date in one phase only")
edges3 <- geom_vline(xintercept = c(ph3$lo, ph3$hi), colour = "grey70",
                     linetype = "dashed", linewidth = 0.3)

strip3 <- ggplot(ph3) +
  edges3 +
  geom_rect(aes(xmin = lo, xmax = hi, ymin = K - phase, ymax = K - phase + 0.8),
            fill = "grey85", colour = "grey30", linewidth = 0.3) +
  geom_text(aes(x = mid, y = K - phase + 0.4, label = paste0("P", phase)),
            size = 3.2) +
  labs(x = NULL, y = NULL) +
  strip_theme

rows3 <- ggplot(ov) +
  edges3 +
  geom_segment(aes(x = Start_date, xend = End_date, y = row, yend = row,
                   colour = "recorded phase window"), linewidth = 2.2) +
  geom_point(aes(True_date, row, shape = where), colour = ACCENT, fill = "white",
             size = 2.2, stroke = 0.9) +
  scale_colour_manual(values = "grey55") +
  scale_shape_manual(values = c("true date in one phase only" = 16,
                                "true date where two phases overlap" = 21)) +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  labs(x = "calendar year", y = ROW_LAB) +
  base_theme

save_fig(strip3, rows3, "Case 3: overlapping phases",
         "A find where two phases overlap is given either phase at random",
         "case3_dating.png", c(0, 1000), strip_h = 1.8)

# Case 4: runs of adjacent abutting phases

mg   <- simulate_merged(N_SHOW + 4, 8, 0.02, 2, K, ALPHA, merge_max = 3,
                        seed = PHASE_SEED)
b4   <- attr(mg, "bounds")
ph4  <- data.frame(phase = seq_len(K), lo = b4[-(K + 1)], hi = b4[-1])
span <- mg$Phases
mg$span <- factor(paste("window of", span, ifelse(span == 1, "phase", "phases")),
                  levels = paste("window of", 1:3, c("phase", "phases", "phases")))
mg   <- sort_rows(mg)
edges4 <- geom_vline(xintercept = b4, colour = "grey70", linetype = "dashed",
                     linewidth = 0.3)

strip4 <- ggplot(ph4) +
  edges4 +
  geom_rect(aes(xmin = lo, xmax = hi, ymin = 0, ymax = 0.8,
                fill = factor(phase %% 2)), colour = "grey30", linewidth = 0.3) +
  geom_text(aes(x = (lo + hi) / 2, y = 0.4, label = paste0("P", phase)), size = 3.2) +
  scale_fill_manual(values = c("grey90", "grey75"), guide = "none") +
  scale_y_continuous(limits = c(-0.3, 1.1)) +
  labs(x = NULL, y = NULL) +
  strip_theme

rows4 <- ggplot(mg) +
  edges4 +
  geom_segment(aes(x = Start_date, xend = End_date, y = row, yend = row,
                   colour = span), linewidth = 2.2) +
  geom_point(aes(True_date, row, shape = TRUE_LAB), fill = ACCENT, colour = "white",
             size = 2.4, stroke = 0.6) +
  scale_colour_manual(values = c("grey20", "grey55", "grey80"), drop = FALSE) +
  scale_shape_manual(values = 21) +
  labs(x = "calendar year", y = ROW_LAB) +
  base_theme

save_fig(strip4, rows4, "Case 4: merged phases",
         "Each find is dated to 1, 2 or 3 adjacent phases",
         "case4_dating.png", c(0, 1000))
