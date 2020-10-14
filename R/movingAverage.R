#' Computes the Moving Average of a Single Time-Series
#'
#' Computes the moving average about a time-series defined by a specified number of points.
#'
#' @param x a \code{vector} of \code{numeric} time-series expression values.
#' @param n a \code{numeric} specifying the number of points to use in the moving average. Default \code{n = 3}.
#' @param centered a \code{logical} scalar. Should the moving average be centered about the current points? Default \code{TRUE} (i.e. average of current point (\code{p}) with  \code{p - n/2} and \code{p + n/2}).
#'
#' @return a \code{vector} containing the smoothed \code{numeric} moving average time-series expression values.
#'
#' @seealso \code{\link{movingAverageDF}}
#'
#' @export
#'
movingAverage <- function(x, n = 3, centered = TRUE) {
  if (centered) {
    before <- floor((n - 1) / 2)
    after <- ceiling((n - 1) / 2)
  } else {
    before <- n - 1
    after <- 0
  }

  # Track the sum and count of number of non-NA items
  s <- rep(0, length(x))
  count <- rep(0, length(x))

  # Add the centered data
  new <- x
  # Add to count list wherever there isn't a
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new

  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new <- c(rep(NA, i), x[1:(length(x) - i)])

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i + 1
  }

  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new <- c(x[(i + 1):length(x)], rep(NA, i))

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i + 1
  }

  # return sum divided by count
  return(s / count)
}
