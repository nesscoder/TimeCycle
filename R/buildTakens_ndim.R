#' Generates Takens' Embedding For a Time-Series
#'
#' Generates the Takens' embedding for a time-series given a specified delay (i.e. lag) and dimension (i.e. 2-D or 3-D).
#'
#' @param x a \code{vector} of \code{numeric} time-series expression values.
#' @param dim a \code{numeric} specifying the dimension to use for in the time-delayed embedding (i.e. 2-D or 3-D).
#' @param delay a \code{numeric} specifying the lag to use for in the n-dimensional time delayed embedding specified by \code{dim}.
#'
#' @seealso \code{\link{getPersistence}}.
#' @return a \code{data.frame} of the n-dimensional Takens' embedding. Columns defined from (1-D to n-D).
#' @export
#'
#'
buildTakens_ndim <- function(x, dim, delay = 1) {
  n <- length(x) - (dim - 1) * delay
  X <- seq_along(x)
  if (n <= 0) {
    stop("Insufficient observations for the requested embedding")
  }
  out <- matrix(rep(X[seq_len(n)], dim), ncol = dim)
  out[, -1] <- out[, -1, drop = FALSE] +
    rep(seq_len(dim - 1) * delay, each = nrow(out))

  return(out)
}
