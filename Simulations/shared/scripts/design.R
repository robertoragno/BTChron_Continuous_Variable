# The part of a design grid that is the same in every case.
#
# Each case stacks its own sweep blocks, then calls add_nuisance() on the
# result. Everything below that point - the ids, the seeds, the per-dataset
# nuisance parameters - was written out three times before, identically.

#' Give every dataset an id, a seed, and its own nuisance parameters.
#'
#' Called once, on the stacked design, so dataset_id is unique across the whole
#' grid; seq_len() per block would have made the blocks collide. The id doubles
#' as the simulation seed.
#'
#' Nuisance parameters vary per dataset, so every cell is a recovery study
#' rather than one configuration repeated. Zero-slope cells are exactly zero,
#' not a small random draw: ratio summaries blow up when divided by a near-zero
#' true slope.
#'
#' Draws three runif() in a fixed order, so a case that draws more of its own
#' afterwards (Case 4 draws K and alpha_conc) keeps a reproducible stream.
add_nuisance <- function(design) {
  design$dataset_id <- seq_len(nrow(design))
  design$seed       <- design$dataset_id

  design$intercept <- runif(nrow(design), 2, 15)
  design$sigma     <- runif(nrow(design), 0.5, 4)
  design$slope     <- ifelse(design$slope_condition == "zero", 0,
                             runif(nrow(design), -0.03, 0.03))
  design
}
