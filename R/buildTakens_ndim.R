#' Generates Taken Embedding
#'
#' Generates the Takens Embedding given a specified delay and dimension
#'
#' @param x int
#' @param dim int
#' @param delay int
#' @param indices boolean
#' @param as.embed boolean
#'
#' @return
#' @export
#'
#'
buildTakens_ndim <- function(x, dim, delay = 1, indices = FALSE, as.embed = FALSE) {
  n <- length(x) - (dim - 1) * delay
  X <- seq_along(x)
  if (n <= 0) {
    stop("Insufficient observations for the requested embedding")
  }
  out <- matrix(rep(X[seq_len(n)], dim), ncol = dim)
  out[, -1] <- out[, -1, drop = FALSE] +
    rep(seq_len(dim - 1) * delay, each = nrow(out))
  if (as.embed) {
    out <- out[, rev(seq_len(ncol(out)))]
  }
  if (!indices) {
    out <- matrix(x[out], ncol = dim)
  }

  return(out)
}
