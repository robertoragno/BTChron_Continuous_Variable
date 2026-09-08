# Case 1 overview figure: what a calibrated date looks like as a weight row, and
# what the models do with a whole dataset of them.
#
# Panel A is one plateau date. Panel B is a 50-date dataset from the same window,
# fitted the same way 00_check_fit.R does, same seed. Each date in Panel B is
# drawn as its whole calibrated distribution, a slab sitting at that date's y
# value, because that whole distribution is what the date-marginalised model
# consumes - not a single summarised year.
#
# Greyscale with the single dark-red accent for the model under test, matching
# the recovery figures (shared/scripts/recovery_summary.R).

library(here)
suppressMessages(library(rcarbon))
library(cmdstanr)
library(ggplot2)
library(patchwork)

source(here("Simulations", "shared", "scripts", "weight_rows.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))
source(here("Simulations", "Sim_Case1_Radiocarbon", "scripts", "simulate.R"))

WINDOW    <- c(-800, -400)
LAB_ERROR <- 25
N         <- 50
STEP      <- 5
INTERCEPT <- 8
SLOPE     <- 0.02
SIGMA     <- 2
OUT <- here("Simulations", "Sim_Case1_Radiocarbon", "figures", "overview.png")

ACCENT   <- "#780000"   # the model under test, as in recovery_summary.R
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

# ---- Panel A: one plateau date as a weight row ----

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

# ---- Panel B: a dataset, and what the models recover ----

pad  <- c14_grid_pad(WINDOW, LAB_ERROR)
grid <- seq(WINDOW[1] - pad, WINDOW[2] + pad, by = STEP)
sim  <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, WINDOW, LAB_ERROR, seed = 1)

x   <- calibrate(sim$CRA, errors = sim$Error, calCurves = "intcal20",
                 calMatrix = TRUE, verbose = FALSE)
calB <- calmatrix_to_calendar(x)
w   <- calibrated_rows(calB, grid)
med_bcad <- 1950 - summary(x)$MedianBP

x_pred <- seq(WINDOW[1], WINDOW[2], length.out = 50)
fit <- cmdstan_model(here("Simulations", "shared", "models",
                          "marginal_date.stan"))$sample(
  data = list(N = N, y = sim$Value, n_years = length(grid), grid_year = grid,
              n_weights = w$n_weights,
              log_year_prob_packed = w$log_year_prob_packed,
              row_first_year = w$row_first_year, row_n_years = w$row_n_years,
              time_ref_min = min(grid), time_ref_range = max(grid) - min(grid),
              N_pred = length(x_pred), x_pred = x_pred),
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, seed = 1, refresh = 0)

mu   <- fit$draws("mu_pred", format = "draws_matrix")
band <- data.frame(x = x_pred,
                   lo = apply(mu, 2, quantile, 0.05),
                   hi = apply(mu, 2, quantile, 0.95))

ols        <- coef(lm(sim$Value ~ med_bcad))
slope_marg <- median(fit$draws("slope_original", format = "draws_matrix")[, 1])
pts        <- data.frame(med = med_bcad, y = sim$Value)

lv <- c("true trend", "date-marginalised", "regression on calibrated medians")
fitlines <- rbind(
  data.frame(kind = lv[1], x = x_pred, y = INTERCEPT + SLOPE * x_pred),
  data.frame(kind = lv[2], x = x_pred, y = apply(mu, 2, median)),
  data.frame(kind = lv[3], x = x_pred, y = ols[1] + ols[2] * x_pred))
fitlines$kind <- factor(fitlines$kind, levels = lv)

pB <- ggplot(fitlines, aes(x, y)) +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE,
              fill = ACCENT, alpha = 0.12) +
  geom_point(data = pts, aes(med, y), inherit.aes = FALSE,
             colour = "grey50", size = 1.3, alpha = 0.75) +
  geom_line(aes(colour = kind, linetype = kind), linewidth = 0.8) +
  scale_colour_manual(values = setNames(c(INK, ACCENT, "grey50"), lv)) +
  scale_linetype_manual(values = setNames(c("solid", "solid", "dashed"), lv)) +
  bce_axis +
  labs(title = "B  Fifty dates from the same window",
       subtitle = sprintf(paste("slope: true %.4f  |  date-marginalised %.4f  |",
                                 "calibrated medians %.4f"),
                          SLOPE, slope_marg, ols[2]),
       caption = paste("Points sit at each date's calibrated median; the",
                       "date-marginalised model reads the whole distribution.",
                       "\nBand: the date-marginalised 90% interval for the trend."),
       x = "calendar year", y = "y") +
  panel_theme

ggsave(OUT, pA / pB, width = 7.4, height = 7.8, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")

# ---- dataset anatomy: plateau vs control, in the shared idiom ----
#
# The flat-window cases get a segment from Start_date to End_date. Radiocarbon
# has no window, so the segment here is each date's 95% HPD envelope - the same
# span the midpoint model collapses the date to - and the grey square is the
# calibrated median. free_x because the two windows are 800 yr apart and would
# not share an axis.

anat_windows <- list("plateau (800-400 BCE)"  = c(-800, -400),
                     "control (1600-1200 BCE)" = c(-1600, -1200))

anat <- do.call(rbind, Map(function(lab, win) {
  s   <- simulate_c14(N, INTERCEPT, SLOPE, SIGMA, win, LAB_ERROR, seed = 1)
  cm  <- calmatrix_to_calendar(
    calibrate(s$CRA, errors = s$Error, calCurves = "intcal20",
              calMatrix = TRUE, verbose = FALSE))
  g   <- seq(win[1] - c14_grid_pad(win, LAB_ERROR),
             win[2] + c14_grid_pad(win, LAB_ERROR), by = STEP)
  smry <- calibrated_summaries(cm, g)
  data.frame(panel = lab, Start_date = smry$start, End_date = smry$end,
             True_date = s$True_date, Value = s$Value,
             Point_date = smry$median)
}, names(anat_windows), anat_windows))
anat$panel <- factor(anat$panel, levels = names(anat_windows))

trend <- data.frame(panel = factor(names(anat_windows),
                                   levels = names(anat_windows)),
                    intercept = INTERCEPT, slope = SLOPE)

ANAT <- here("Simulations", "Sim_Case1_Radiocarbon", "figures",
             "dataset_anatomy.png")
ggsave(ANAT,
       dataset_anatomy(anat, trend = trend, free_x = TRUE,
                       point_label = "calibrated median",
                       span_label  = "95% HPD envelope",
                       title = "Case 1: one dataset per window"),
       width = 9, height = 3.8, dpi = 300, bg = "white")
cat("wrote", ANAT, "\n")
