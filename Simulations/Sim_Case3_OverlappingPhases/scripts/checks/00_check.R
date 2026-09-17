# What the Case 3 generator produces, before anything is fitted to it.
# Prints the midpoint slope ratio (plain least squares, no Stan) and writes
# figures/checks/check_depletion.png.

library(here)
library(ggplot2)

source(here("Simulations", "Sim_Case3_OverlappingPhases", "scripts", "simulate.R"))

OVERLAPS  <- c(0, 0.25, 0.5)
ASSIGN_PS <- c(0.5, 0.35, 0.2)
SLOPE     <- 0.02
N_REP     <- 200   # datasets per setting

midpoint <- function(d) (d$Start_date + d$End_date) / 2

# Midpoint slope ratio: 1 means the midpoint recovers the slope
settings <- expand.grid(overlap = OVERLAPS, assign_p = ASSIGN_PS)
settings <- settings[settings$overlap > 0 | settings$assign_p == 0.5, ]
settings$slope_ratio <- sapply(seq_len(nrow(settings)), function(g) {
  ratios <- sapply(seq_len(N_REP), function(i) {
    d <- simulate_overlap(200, 8, SLOPE, 1, K = sample(3:10, 1),
                          alpha_conc = 10^(-1 + 2 * rbeta(1, 2, 3.5)),
                          overlap = settings$overlap[g],
                          assign_p = settings$assign_p[g], seed = 1000 + i)
    coef(lm(d$Value ~ midpoint(d)))[2] / SLOPE
  })
  # NA when every find lands in one phase, so the midpoint does not vary
  median(ratios, na.rm = TRUE)
})
print(settings, row.names = FALSE, digits = 3)

# True dates kept by one wide phase at overlap 0.5, for each assign_p
K <- 5
bounds <- attr(simulate_overlap(10, 8, SLOPE, 1, K, alpha_conc = 2, overlap = 0,
                                seed = 7), "bounds")
widest <- which.max(diff(bounds)[2:(K - 1)]) + 1   # widest phase not at an end

kept <- do.call(rbind, lapply(ASSIGN_PS, function(ap) {
  d <- simulate_overlap(40000, 8, SLOPE, 1, K, alpha_conc = 2, overlap = 0.5,
                        assign_p = ap, seed = 7)
  phase <- findInterval(midpoint(d), bounds, rightmost.closed = TRUE,
                        all.inside = TRUE)
  d <- d[phase == widest, ]
  data.frame(assign_p = paste("assign_p", ap), True_date = d$True_date,
             midpoint = midpoint(d)[1], mean_date = mean(d$True_date))
}))

lines <- unique(kept[, c("assign_p", "midpoint", "mean_date")])
lines <- rbind(data.frame(assign_p = lines$assign_p, x = lines$midpoint,
                          line = "window midpoint"),
               data.frame(assign_p = lines$assign_p, x = lines$mean_date,
                          line = "mean true date"))

p <- ggplot(kept, aes(True_date)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "grey30", linewidth = 0.2) +
  geom_vline(data = lines, aes(xintercept = x, linetype = line, colour = line),
             linewidth = 0.6) +
  scale_colour_manual(values = c("window midpoint" = "grey30",
                                 "mean true date" = "#780000")) +
  scale_linetype_manual(values = c("window midpoint" = "dashed",
                                   "mean true date" = "solid")) +
  facet_wrap(~ assign_p, ncol = 1) +
  labs(x = "True date (yr)", y = "Count", colour = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), legend.position = "top")

out <- here("Simulations", "Sim_Case3_OverlappingPhases", "figures", "checks",
            "check_depletion.png")
ggsave(out, p, width = 7.5, height = 6.5, dpi = 300, bg = "white")
cat("wrote", out, "\n")
