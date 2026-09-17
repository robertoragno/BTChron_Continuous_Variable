# What the Case 2 generator produces, before anything is fitted to it.
# Writes figures/checks/check.png and figures/dataset_anatomy.png.

library(here)
library(ggplot2)
library(patchwork)

source(here("Simulations", "Sim_Case2_Typochronology", "scripts", "simulate.R"))
source(here("Simulations", "shared", "scripts", "recovery_summary.R"))

N_BIG <- 20000
fig   <- function(name) here("Simulations", "Sim_Case2_Typochronology", "figures", name)
theme_check <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"))

sim <- function(n, ..., seed = 1) simulate_typo(n, 8, 0.02, 1, ..., seed = seed)

# Every window must contain its own true date
d <- sim(N_BIG)
stopifnot(all(d$True_date >= d$Start_date & d$True_date <= d$End_date))

# 1. True dates: even, and denser at the end with growth_ratio = 4
dates <- rbind(data.frame(deposition = "even", t = sim(N_BIG)$True_date),
               data.frame(deposition = "growth ratio 4",
                          t = sim(N_BIG, growth_ratio = 4)$True_date))
p_dates <- ggplot(dates, aes(t)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", linewidth = 0.2) +
  facet_wrap(~ deposition) +
  labs(title = "True dates", x = "Year", y = "Count") + theme_check

# 2. Share of coarsely dated finds along the period, with and without the trend
share <- rbind(data.frame(trend = "no trend", d = I(list(sim(N_BIG)))),
               data.frame(trend = "precision trend 0.6",
                          d = I(list(sim(N_BIG, precision_trend = 0.6)))))
share <- do.call(rbind, lapply(seq_len(nrow(share)), function(i)
  data.frame(trend = share$trend[i], t = share$d[[i]]$True_date,
             coarse = share$d[[i]]$Coarse)))
p_share <- ggplot(share, aes(t, as.numeric(coarse))) +
  stat_summary_bin(bins = 16, fun = mean, geom = "point") +
  facet_wrap(~ trend) + ylim(0, 1) +
  labs(title = "Share of coarsely dated finds along the period",
       x = "Year", y = "Share coarsely dated") + theme_check

# 3. Where the true date sits inside its window: must be flat
p_pos <- ggplot(data.frame(pos = (d$True_date - d$Start_date) / (d$End_date - d$Start_date)),
                aes(pos)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "black", linewidth = 0.2) +
  labs(title = "True date within its window (must be flat)",
       x = "(true - start) / width", y = "Count") + theme_check

ggsave(fig(file.path("checks", "check.png")), p_dates / p_share / p_pos,
       width = 9, height = 9, dpi = 300, bg = "white")

# One example dataset per proportion of coarsely dated finds
shares <- c(0, 0.5, 1)
labels <- paste0(100 * shares, "% coarsely dated")
anat <- do.call(rbind, lapply(seq_along(shares), function(i) {
  s <- sim(60, prop_coarse_samples = shares[i], seed = 7)
  data.frame(panel = factor(labels[i], levels = labels),
             Start_date = s$Start_date, End_date = s$End_date,
             True_date = s$True_date, Value = s$Value,
             Point_date = (s$Start_date + s$End_date) / 2)
}))
trend <- data.frame(panel = factor(labels, levels = labels),
                    intercept = 8, slope = 0.02)
ggsave(fig("dataset_anatomy.png"),
       dataset_anatomy(anat, trend = trend,
                       title = "Case 2: one dataset per proportion of coarsely dated finds"),
       width = 9, height = 3.6, dpi = 300, bg = "white")

cat("checks passed, figures written\n")
