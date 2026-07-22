# Real-data linear case study: prepare the GINI subset and the exploratory panel.
#
# Subset: Great Britain, nuclear-family houses, positive roofed area, with a
# dated window (EndHouse > BeginHouse) and midpoint in [-100, 1000].
# Response is log(RoofedArea); dating window is [BeginHouse, EndHouse].
#
# Produces: gini_gb_filtered.csv, dataset_1_exploratory_panel.png

library(here)
library(tidyverse)
library(patchwork)

data_dir   <- here("Real_Data", "dataset_1", "data")
raw_path   <- file.path(data_dir, "gini_database_all_records_20240721.csv")
clean_path <- file.path(data_dir, "gini_gb_filtered.csv")
figure_path <- function(name) here("Real_Data", "dataset_1", "figures", name)

panel_theme <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.margin = margin(5, 10, 2, 10))

# 1. Filter ---------------------------------------------------------------------

raw <- read_csv(raw_path, show_col_types = FALSE)

gb_data <- raw %>%
  filter(Region == "Great Britain", RoofedArea > 0,
         Housetypedesc == "Nuclear family", EndHouse > BeginHouse) %>%
  mutate(Midpoint = (BeginHouse + EndHouse) / 2) %>%
  filter(Midpoint >= -100, Midpoint <= 1000) %>%
  transmute(ID = row_number(),
            Site,
            Start_date = BeginHouse,
            End_date   = EndHouse,
            Midpoint,
            RoofedArea,
            Value = log(RoofedArea))

write_csv(gb_data, clean_path)
cat("Filtered subset:", nrow(gb_data), "houses\n")

# 2. Exploratory panel ----------------------------------------------------------
# Same three-panel layout as the simulation, without a ground-truth line:
# real data has no known trend.

year_breaks <- seq(-200, 1000, 200)

exploratory <- gb_data %>%
  arrange(Midpoint) %>%
  mutate(row_position = row_number())

values_histogram <- ggplot(exploratory, aes(Value)) +
  geom_histogram(binwidth = 0.3, fill = "grey70", colour = "black", linewidth = 0.3) +
  labs(title = expression(bold("A.") ~ "Distribution of values"),
       x = "log(Roofed area)", y = "Count") + panel_theme

date_ranges_plot <- ggplot(exploratory, aes(y = row_position)) +
  geom_segment(aes(x = Start_date, xend = End_date, yend = row_position),
               linewidth = 0.1, colour = "grey40", alpha = 0.5) +
  geom_point(aes(x = Midpoint), shape = 15, size = 0.2, colour = "black") +
  scale_x_continuous(breaks = year_breaks, expand = c(0.01, 0)) +
  labs(title = expression(bold("B.") ~ "Recorded dating ranges"), x = "Year",
       y = "House (ordered)") +
  panel_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

values_over_time <- ggplot(exploratory, aes(Midpoint, Value)) +
  geom_linerange(aes(xmin = Start_date, xmax = End_date), linewidth = 0.1,
                 colour = "grey60", alpha = 0.3) +
  geom_point(shape = 15, size = 0.4, colour = "black", alpha = 0.4) +
  scale_x_continuous(breaks = year_breaks, expand = c(0.01, 0)) +
  labs(title = expression(bold("C.") ~ "Values across time"), x = "Year",
       y = "log(Roofed area)") + panel_theme

ggsave(figure_path("dataset_1_exploratory_panel.png"),
       (values_histogram | date_ranges_plot | values_over_time) +
         plot_layout(widths = c(0.8, 1.8, 0.8)),
       width = 14, height = 4.5, dpi = 300, bg = "white")

cat("Wrote dataset_1_exploratory_panel.png\n")
