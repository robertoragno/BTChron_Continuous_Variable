#!/usr/bin/env Rscript
setwd("/home/rr673/R_projects/BTChron_Paper_1")

cat("══ Fitting model ══════════════════════════════════════════════════════\n")
source("Simulations/Sim_Linear/scripts/03_fit_linear_model.R")

cat("\n══ Generating result plots ═════════════════════════════════════════════\n")
source("Simulations/Sim_Linear/scripts/04_plot_results.R")

cat("\n══ Done ════════════════════════════════════════════════════════════════\n")
