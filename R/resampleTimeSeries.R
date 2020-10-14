#' Generates a \code{data.frame} of Resampled Time-Series
#'
#' Converts a \code{data.frame} of \code{numeric} values into a single vector and generates a random resampling with dimension
#' \code{numRsmps} by \code{numTP}.
#'
#' @param distTP a \code{data.frame} of \code{numeric} values.
#' @param numTP a \code{numeric} specifying the number of columns in the outputted \code{data.frame}.
#' @param numRsmps a \code{numeric} specifying the number of rows in the outputted \code{data.frame}.
#'
#' @return a \code{data.frame} of randomly sampled \code{numerics} time-series (row = \code{numRsamps} \emph{x} col = \code{numTP}).
#'
#' @seealso \code{\link{nullResampling}}
#'
#' @export
#'
resampleTimeSeries <- function(distTP, numTP, numRsmps) {
  dist <- as.numeric(as.vector(unlist(distTP)))
  return(as.data.frame(matrix(data = sample(dist, size = numTP * numRsmps, replace = T), ncol = numTP, nrow = numRsmps)))
}
