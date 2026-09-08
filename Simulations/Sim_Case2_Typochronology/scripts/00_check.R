# What the Case 2 generator produces, before anything is fitted to it.
# Meant to be stepped through: edit the settings below, re-source, look at the
# figure. Running the whole file with Rscript does the same thing.

library(here)
library(ggplot2)
library(patchwork)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

# ---- settings to play with ----
PERIOD_START <- 100     # matches Sim_Linear's TMIN/TMAX, so Case 3 is comparable
PERIOD_END   <- 900
MIN_FRAC     <- 0.05    # narrowest window, as a share of the period
MAX_FRAC     <- 0.25    # widest window; beyond ~1/3 the trend starts flattening
GRID         <- 25      # window boundaries are multiples of this
STRENGTH     <- 4       # how lopsided the fine/coarse mixes are (1 = flat)
SKEW_STRONG  <- 4       # Beta shape for the 2b panels
N_BIG        <- 20000   # for the distribution panels
# -------------------------------

widths <- typo_widths(PERIOD_START, PERIOD_END, MIN_FRAC, MAX_FRAC, GRID)
cat("width set:", paste(widths, collapse = ", "), "\n")
cat("as % of period:",
    paste0(round(100 * widths / (PERIOD_END - PERIOD_START)), "%", collapse = ", "), "\n")

figure_path <- function(name) here("Simulations", "Sim_Case2_Typochronology", "figures", name)
panel_theme  <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "top", legend.title = element_blank())

sim <- function(n, mix = "even", skew = 1, seed = 1) {
  simulate_typo(n, intercept = 8, slope = 0.02, sigma = 1,
                period_start = PERIOD_START, period_end = PERIOD_END,
                widths = widths,
                width_prob = typo_width_prob(widths, mix, STRENGTH),
                skew_shape = skew, grid = GRID, seed = seed)
}

# The generator must never hand out a window that fails to contain its own date,
# or boundaries off the grid. Stop rather than plot if either breaks.
for (m in c("even", "fine", "coarse")) {
  d <- sim(5000, mix = m, skew = SKEW_STRONG)
  stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date),
            all(d$Start_date %% GRID == 0), all(d$End_date %% GRID == 0))
}
cat("containment and grid checks passed for all three mixes\n")

# 1. Dating resolution: what fine / even / coarse actually look like.
mix_df <- do.call(rbind, lapply(c("fine", "even", "coarse"), function(m) {
  d <- sim(N_BIG, mix = m)
  data.frame(mix = m, Width = d$End_date - d$Start_date)
}))
mix_df$mix <- factor(mix_df$mix, levels = c("fine", "even", "coarse"))

mix_plot <- ggplot(mix_df, aes(factor(Width))) +
  geom_bar(fill = "grey70", colour = "black", linewidth = 0.2) +
  facet_wrap(~ mix) +
  labs(title = "Dating resolution: the fine / coarse factor",
       x = "Window width (yr)", y = "Count") + panel_theme

# 2. Where the true date sits in its window. Uniform for 2a, leaning late for 2b.
frac_df <- rbind(
  data.frame(case = "2a (uniform)",
             frac = with(sim(N_BIG), (True_date - Start_date) / (End_date - Start_date))),
  data.frame(case = sprintf("2b (skew = %d)", SKEW_STRONG),
             frac = with(sim(N_BIG, skew = SKEW_STRONG),
                         (True_date - Start_date) / (End_date - Start_date))))

frac_plot <- ggplot(frac_df, aes(frac)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "black", linewidth = 0.2) +
  facet_wrap(~ case) +
  labs(title = "True date within its window",
       x = "(true - start) / width", y = "Count") + panel_theme

# 3. The two that catch window placement leaking into the x-design. Both flat.
d_even <- sim(N_BIG)
width_by_pos <- ggplot(data.frame(t = d_even$True_date,
                                  w = d_even$End_date - d_even$Start_date),
                       aes(t, w)) +
  stat_summary_bin(bins = 20, fun = mean, geom = "point", size = 1.5) +
  ylim(0, NA) +
  labs(title = "Mean width by position (must be flat)",
       x = "Year", y = "Mean width (yr)") + panel_theme

date_plot <- ggplot(d_even, aes(True_date)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", linewidth = 0.2) +
  ylim(0, NA) +
  labs(title = "True dates (must be flat)", x = "Year", y = "Count") + panel_theme

ggsave(figure_path("check.png"),
       mix_plot / frac_plot / (width_by_pos | date_plot),
       width = 9, height = 9.5, dpi = 300, bg = "white")

# 4. An example dataset small enough to read as a finds table.
example <- sim(12, seed = 4)
example <- example[order(example$Start_date), ]
example$window <- paste0(example$Start_date, "-", example$End_date, " CE")
example$Width  <- example$End_date - example$Start_date
print(example[, c("window", "Width", "True_date", "Value")], row.names = FALSE)

# 4b. The same view the archived Sim_Linear figures gave: one example dataset per
# resolution, each find's window as a segment, its true date against the trend.
# The point-date here is the window midpoint - what the midpoint model sees.
anat <- do.call(rbind, lapply(c("fine", "even", "coarse"), function(m) {
  d <- sim(60, mix = m, seed = 7)
  data.frame(panel = factor(m, levels = c("fine", "even", "coarse")),
             Start_date = d$Start_date, End_date = d$End_date,
             True_date = d$True_date, Value = d$Value,
             Point_date = (d$Start_date + d$End_date) / 2)
}))
trend <- data.frame(panel = factor(c("fine", "even", "coarse"),
                                   levels = c("fine", "even", "coarse")),
                    intercept = 8, slope = 0.02)
ggsave(figure_path("dataset_anatomy.png"),
       dataset_anatomy(anat, trend = trend,
                       title = "Case 2: one dataset per dating resolution"),
       width = 9, height = 3.6, dpi = 300, bg = "white")

# 5. What the midpoint costs, by mix. OLS only, so this is the dating structure
# and not any Stan model. Ratio < 1 means the trend comes out flattened.
for (m in c("fine", "even", "coarse")) {
  r <- sapply(1:200, function(i) {
    s <- sim(200, mix = m, seed = 500 + i)
    coef(lm(s$Value ~ I((s$Start_date + s$End_date) / 2)))[2] / 0.02
  })
  cat(sprintf("%-7s mean width %3.0f yr   midpoint slope ratio %.3f\n",
              m, mean(sim(5000, mix = m)$End_date - sim(5000, mix = m)$Start_date),
              median(r)))
}
