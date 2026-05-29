library(here)
library(tidyverse)
library(patchwork)

sim_data     <- read_csv(here("Simulations", "Sim_GP", "data", "simulated_data.csv"),
                         show_col_types = FALSE)
ground_truth <- read_csv(here("Simulations", "Sim_GP", "data", "ground_truth.csv"),
                         show_col_types = FALSE)

sim_data <- sim_data %>%
  mutate(Midpoint = (Start_date + End_date) / 2) %>%
  arrange(Midpoint) %>%
  mutate(rank = row_number())

# ── Shared theme ─────────────────────────────────────────────────────────────

theme_panel <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin  = margin(5, 10, 2, 10),
    axis.title.x = element_blank()
  )

# ── Panel A: Histogram of values ────────────────────────────────────────────

p_hist <- ggplot(sim_data, aes(x = Value)) +
  geom_histogram(binwidth = 2, fill = "grey70", colour = "black",
                 linewidth = 0.3) +
  labs(
    title = expression(bold("A.") ~ "Distribution of values"),
    y = "Count"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

# ── Panel B: Duration bars sorted by midpoint ───────────────────────────────

sim_data <- sim_data %>%
  mutate(rank_spaced = rank * 1.8)

tick_h <- 0.5

p_dur <- ggplot(sim_data, aes(y = rank_spaced)) +
  geom_segment(aes(x = Start_date, xend = End_date, yend = rank_spaced),
               linewidth = 0.2, colour = "grey40") +
  geom_segment(aes(x = Start_date, xend = Start_date,
                   y = rank_spaced - tick_h, yend = rank_spaced + tick_h),
               linewidth = 0.3, colour = "grey30") +
  geom_segment(aes(x = End_date, xend = End_date,
                   y = rank_spaced - tick_h, yend = rank_spaced + tick_h),
               linewidth = 0.3, colour = "grey30") +
  geom_point(aes(x = Midpoint), shape = 15, size = 0.35, colour = "black") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  labs(
    title = expression(bold("B.") ~ "Sampled date ranges"),
    x = "Year CE",
    y = "Sample (ordered)"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# ── Panel C: Value vs time with true trend ─────────────────────────────────

p_val <- ggplot(sim_data, aes(x = Midpoint, y = Value)) +
  geom_linerange(aes(xmin = Start_date, xmax = End_date),
                 linewidth = 0.15, colour = "grey60", alpha = 0.5) +
  geom_point(shape = 15, size = 0.6, colour = "black", alpha = 0.7) +
  geom_line(data = ground_truth, aes(x = Year, y = True_value),
            linewidth = 0.7, colour = "black", linetype = "dashed") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("C.") ~ "Distribution of values across time"),
    x = "Year CE",
    y = "Value"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

# ── Combine ──────────────────────────────────────────────────────────────────

p <- (p_hist / p_dur / p_val) +
  plot_layout(heights = c(0.8, 1.8, 0.8)) +
  plot_annotation(
    caption = "Dashed line (C): true Gaussian bell-curve trend used to generate the simulated values.",
    theme = theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey40"))
  )

ggsave(here("Simulations", "Sim_GP", "figures", "exploratory_panel.png"), p,
       width = 6, height = 9, dpi = 300, bg = "white")

cat("Saved to", here("Simulations", "Sim_GP", "figures", "exploratory_panel.png"), "\n")
