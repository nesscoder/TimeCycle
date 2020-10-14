#' Computes the Fast Fourier Transform of a Time-Series
#'
#' Takes a time-series and sample frequency and computes the fast Fourier transform of the time-series.
#'
#' @param y a \code{vector} of numeric time-series expression values.
#' @param sampFreq a \code{vector} of numeric frequencies.
#'
#' @return a \code{data.frame} harmonic fits by frequency.
#'
#' @seealso \code{\link[stats]{fft}}
#' @export
#'
getFFT <- function(y, sampFreq) {
  N <- length(y)
  fk <- stats::fft(y) / N # normalize Data
  fk <- 2 * fk[1:((length(fk) / 2) + 1)] # DC comp + half of positives
  freq <- (0:(N - 1)) * sampFreq / N
  freq <- freq[(1:(length(fk)))]
  return(data.frame(fur = fk, freq = freq))
}
