# Case 1 illustration figures.
#   figures/checks/overview.png  A: one calibrated date on the Hallstatt plateau
#                                B: a 50-date dataset, full distribution vs
#                                   regression on calibrated medians
#   figures/dataset_anatomy.png  one dataset per window, before any fitting

library(here)
suppressMessages(library(rcarbon))
library(cmdstanr)
library(ggplot2)
library(patchwork)

source(here("Simulations", "shared", "dating", "case1_radiocarbon.R"))

WINDOW    <- c(-800, -400)
LAB_ERROR <- 25
N         <- 50
STEP      <- 5
INTERCEPT <- 8
SLOPE     <- 0.02
SIGMA     <- 2
OUT <- here("Simulations", "time_as_predictor", "Case1_Radiocarbon", "figures", "checks", "overview.png")

ACCENT   <- "#780000"   # the full distribution, as in the recovery figures
INK      <- "grey15"
INK_SOFT <- "grey35"

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = INK_SOFT),
        plot.caption = element_text(size = 8, colour = INK_SOFT, hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        legend.key.width = unit(1.4, "lines"))

bce_axis <- scale_x_continuous(
  labels = function(x) ifelse(x < 0, paste0(abs(x), " BCE"), x),
  expand = expansion(mult = 0.03))

# Panel A: one calibrated date on the plateau

demo <- calibrate(2450, errors = 25, calCurves = "intcal20",
                  calMatrix = TRUE, verbose = FALSE)
cal  <- calmatrix_to_calendar(demo)
keep <- cal$years >= -820 & cal$years <= -340
dA   <- data.frame(year = cal$years[keep], prob = cal$prob[keep, 1])

h   <- hpdi(demo, credMass = 0.95)[[1]]
env <- c(1950 - max(h[, "startCalBP"]), 1950 - min(h[, "endCalBP"]))
med <- 1950 - summary(demo)$MedianBP

env_df  <- data.frame(xmin = env[1], xmax = env[2],
                      fill = "95% HPD envelope (4 disjoint chunks)")
marks_A <- data.frame(
  x    = c(med, mean(env)),
  what = factor(c("calibrated median", "HPD-envelope midpoint"),
                levels = c("calibrated median", "HPD-envelope midpoint")))

pA <- ggplot(dA, aes(year, prob)) +
  geom_rect(data = env_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill)) +
  geom_area(fill = "grey72", colour = "grey45", linewidth = 0.3) +
  geom_vline(data = marks_A, aes(xintercept = x, colour = what, linetype = what),
             linewidth = 0.7) +
  scale_fill_manual(values = c("95% HPD envelope (4 disjoint chunks)" = "grey86")) +
  scale_colour_manual(values = c("calibrated median" = ACCENT,
                                 "HPD-envelope midpoint" = "grey30")) +
  scale_linetype_manual(values = c("calibrated median" = "solid",
                                   "HPD-envelope midpoint" = "dashed")) +
  guides(fill = guide_legend(order = 1),
         colour = guide_legend(order = 2), linetype = guide_legend(order = 2)) +
  bce_axis +
  labs(title = "A  One radiocarbon date on the Hallstatt plateau",
       subtitle = "2450 +/- 25 BP calibrated against IntCal20",
       x = NULL, y = "probability per year") +
  panel_theme

# Panel B: a dataset, fitted with the full distribution as in 03_recovery_study.R

pad  <- c14_grid_pad(WINDOW, LAB_ERROR)
grid <- seq(WINDOW[1] - pad, WINDOW[2] + pad, by = STEP)
sim  <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, WINDOW, LAB_ERROR, seed = 1)

# Calibrate over the whole grid, then add up the yearly probabilities into
# 5-year cells: one row per date, one column per grid year
cal <- calibrate(sim$CRA, sim$Error, calMatrix = TRUE, verbose = FALSE,
                 timeRange = 1950 - c(min(grid) - STEP, max(grid) + STEP))
years <- 1950 - as.numeric(rownames(cal$calmatrix))
cell <- round((years - grid[1]) / STEP) + 1
on.grid <- cell >= 1 & cell <= length(grid)
prob <- t(rowsum(cal$calmatrix[on.grid, ], cell[on.grid]))
prob <- prob / rowSums(prob)
first <- rep(NA, N)
last <- rep(NA, N)
med_bcad <- rep(NA, N)
for (i in 1:N)
{
  cum <- cumsum(prob[i, ])
  first[i] <- which(cum >= 0.00005)[1]
  last[i] <- which(cum >= 0.99995)[1]
  med_bcad[i] <- grid[which(cum >= 0.5)[1]]
}

fit <- cmdstan_model(here("Simulations", "shared", "models", "linear_dates.stan"))$sample(
  data = list(n = N, n_years = length(grid), y = sim$Value,
              year = matrix(grid, nrow = N, ncol = length(grid), byrow = TRUE),
              log_p = log(prob), first = first, last = last, centre = mean(WINDOW)),
  chains = 4, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500,
  adapt_delta = 0.95, seed = 1, refresh = 0, show_messages = FALSE)
draws <- as.data.frame(fit$draws(c("intercept", "slope"), format = "draws_df"))

x_pred <- seq(WINDOW[1], WINDOW[2], length.out = 50)
mu <- sapply(x_pred, function(x) draws$intercept + draws$slope * x)
band <- data.frame(x = x_pred, lo = apply(mu, 2, quantile, 0.05),
                   hi = apply(mu, 2, quantile, 0.95), what = "full distribution, 90% interval")

ols        <- coef(lm(sim$Value ~ med_bcad))
slope_marg <- median(draws$slope)
pts        <- data.frame(med = med_bcad, y = sim$Value, what = "calibrated median of each date")

lv <- c("true trend", "full distribution", "regression on calibrated medians")
fitlines <- rbind(
  data.frame(kind = lv[1], x = x_pred, y = INTERCEPT + SLOPE * x_pred),
  data.frame(kind = lv[2], x = x_pred, y = apply(mu, 2, median)),
  data.frame(kind = lv[3], x = x_pred, y = ols[1] + ols[2] * x_pred))
fitlines$kind <- factor(fitlines$kind, levels = lv)

pB <- ggplot(fitlines, aes(x, y)) +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi, fill = what), inherit.aes = FALSE,
              alpha = 0.12) +
  geom_point(data = pts, aes(med, y, shape = what), inherit.aes = FALSE,
             colour = "grey50", size = 1.3, alpha = 0.75) +
  geom_line(aes(colour = kind, linetype = kind), linewidth = 0.8) +
  scale_colour_manual(values = setNames(c(INK, ACCENT, "grey50"), lv)) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), lv)) +
  scale_fill_manual(values = c("full distribution, 90% interval" = ACCENT)) +
  scale_shape_manual(values = c("calibrated median of each date" = 16)) +
  guides(colour = guide_legend(order = 1, nrow = 1), linetype = guide_legend(order = 1, nrow = 1),
         fill = guide_legend(order = 2, theme = theme(legend.key.height = unit(0.8, "lines"))),
         shape = guide_legend(order = 2)) +
  bce_axis +
  labs(title = "B  Fifty dates from the same window",
       subtitle = sprintf(paste("slope: true %.4f  |  full distribution %.4f  |",
                                 "calibrated medians %.4f"),
                          SLOPE, slope_marg, ols[2]),
       x = "calendar year", y = "y") +
  panel_theme +
  theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"), legend.margin = margin(0, 0, 0, 0))

ggsave(OUT, pA / pB, width = 7.4, height = 7.8, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")

# Dataset anatomy: one dataset per window, before any fitting. Each date is drawn
# at its value, with its 95% HPD range, its calibrated median and its true date.
# The two windows get their own x axis.

anat_windows <- list("Hallstatt plateau (800–400 BCE)" = c(-800, -400),
                     "Steep section (1600–1200 BCE)"   = c(-1600, -1200))
N_ANAT <- 25

anat <- data.frame()
trend <- data.frame()
for (lab in names(anat_windows))
{
  win <- anat_windows[[lab]]
  # Shift the intercept per window so both panels sit at the same height
  icpt <- INTERCEPT - SLOPE * mean(win)
  s <- simulate_c14(N_ANAT, icpt, SLOPE, SIGMA, win, LAB_ERROR, seed = 1)
  g <- seq(win[1] - c14_grid_pad(win, LAB_ERROR), win[2] + c14_grid_pad(win, LAB_ERROR), by = STEP)
  cal <- calibrate(s$CRA, s$Error, calMatrix = TRUE, verbose = FALSE,
                   timeRange = 1950 - c(min(g) - STEP, max(g) + STEP))
  years <- 1950 - as.numeric(rownames(cal$calmatrix))
  cell <- round((years - g[1]) / STEP) + 1
  on.grid <- cell >= 1 & cell <= length(g)
  p <- t(rowsum(cal$calmatrix[on.grid, ], cell[on.grid]))
  p <- p / rowSums(p)
  for (i in 1:N_ANAT)
  {
    # 95% HPD range: the most probable years until they hold 95% of the probability
    cut <- sort(p[i, ], decreasing = TRUE)[which(cumsum(sort(p[i, ], decreasing = TRUE)) >= 0.95)[1]]
    anat <- rbind(anat, data.frame(panel = lab, start = min(g[p[i, ] >= cut]), end = max(g[p[i, ] >= cut]),
                                   median = g[which(cumsum(p[i, ]) >= 0.5)[1]],
                                   true = s$True_date[i], value = s$Value[i]))
  }
  trend <- rbind(trend, data.frame(panel = lab, intercept = icpt, slope = SLOPE))
}
anat$panel <- factor(anat$panel, levels = names(anat_windows))
trend$panel <- factor(trend$panel, levels = names(anat_windows))

p_anat <- ggplot(anat) +
  geom_abline(data = trend, aes(intercept = intercept, slope = slope, linetype = "true trend")) +
  geom_segment(aes(x = start, xend = end, y = value, yend = value, colour = "95% HPD range")) +
  geom_point(aes(x = median, y = value, colour = "calibrated median"), shape = 15, size = 1.2) +
  geom_point(aes(x = true, y = value, colour = "true date"), size = 1) +
  scale_colour_manual(values = c("95% HPD range" = "grey70", "calibrated median" = "grey40",
                                 "true date" = ACCENT)) +
  scale_linetype_manual(values = c("true trend" = "dashed")) +
  facet_wrap(~panel, scales = "free_x") +
  bce_axis +
  labs(title = "Case 1: one dataset per window", x = "calendar year", y = "measured value") +
  panel_theme

ANAT <- here("Simulations", "time_as_predictor", "Case1_Radiocarbon", "figures",
             "dataset_anatomy.png")
ggsave(ANAT, p_anat, width = 9, height = 3.8, dpi = 300, bg = "white")
cat("wrote", ANAT, "\n")
