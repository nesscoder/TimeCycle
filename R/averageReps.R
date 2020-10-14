#' Handles Sample Replicates in Single Time-Series
#'
#' Averages expression values for a single time-series across replicate groups.
#'
#' @param geneExpr a \code{vector} of \code{numeric} time-series expression values with replicates.
#' @param Reps a \code{vector} defining the number of replicates at each time point.
#'
#' @return a \code{vector} of average  \code{numeric} time-series expression values by replicate time-points.
#' @seealso \code{\link{getRepAvgedDataFrame}}
#'
#' @export
#'
#' @examples
#'
#' geneExpr <- c(1, 5, 3, 4, 5, 8)
#' reps <- c(2, 3, 1)
#' averageReps(geneExpr = geneExpr, Reps = reps)
averageReps <- function(geneExpr, Reps) {
  splitBy <- unlist(sapply(1:length(Reps), function(seq) {
    rep(seq, Reps[seq])
  }))
  as.vector(unlist(lapply(split(geneExpr, splitBy), mean)))
}
