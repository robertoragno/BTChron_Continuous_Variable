#' Give every dataset in a design an id, a seed, and its own true trend.
#'
#' Called once on the whole design, so ids are unique; the id is also the seed.
#' Intercept, sigma and slope are drawn per dataset (in this order, which keeps
#' later random draws reproducible). Zero-slope datasets get a slope of exactly 0.
add_true_trend <- function(design) {
  design$dataset_id <- seq_len(nrow(design))
  design$seed       <- design$dataset_id

  design$intercept <- runif(nrow(design), 2, 15)
  design$sigma     <- runif(nrow(design), 0.5, 4)
  design$slope     <- ifelse(design$slope_condition == "zero", 0,
                             runif(nrow(design), -0.03, 0.03))
  design
}
