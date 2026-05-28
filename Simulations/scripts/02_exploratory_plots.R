library(here)
library(tidyverse)
library(patchwork)

sim_data <- read_csv(here("Simulations", "data", "simulated_data.csv"),
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

p_dur <- ggplot(sim_data, aes(y = rank)) +
  geom_segment(aes(x = Start_date, xend = End_date, yend = rank),
               linewidth = 0.15, colour = "grey30") +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("B.") ~ "Sample date ranges (sorted by midpoint)"),
    x = "Year CE",
    y = "Sample (ordered)"
  ) +
  theme_panel +
  theme(
    axis.title.x = element_text(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# ── Panel C: Value on Y axis, date ranges on X ─────────────────────────────

p_val <- ggplot(sim_data) +
  geom_segment(aes(x = Start_date, xend = End_date, y = Value, yend = Value),
               linewidth = 0.15, colour = "grey50", alpha = 0.4) +
  geom_point(aes(x = Midpoint, y = Value), shape = 16, colour = "black",
             size = 0.6, alpha = 0.5) +
  scale_x_continuous(breaks = seq(100, 900, 100), expand = c(0.01, 0)) +
  labs(
    title = expression(bold("C.") ~ "Value vs. date range"),
    x = "Year CE",
    y = "Value"
  ) +
  theme_panel +
  theme(axis.title.x = element_text())

# ── Combine ──────────────────────────────────────────────────────────────────

p <- (p_hist / p_dur / p_val) +
  plot_layout(heights = c(1, 1.4, 1))

out_dir <- here("Simulations", "figures")
ggsave(file.path(out_dir, "exploratory_panel.png"), p,
       width = 6, height = 8, dpi = 300, bg = "white")

cat("Saved to", file.path(out_dir, "exploratory_panel.png"), "\n")
