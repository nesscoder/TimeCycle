#' getFfT
#'
#' computes the fast Fourier Transform of a Time series
#'
#' @param y double
#' @param sampFreq double
#'
#' @return
#' @export
#'
getFFT <- function(y, sampFreq) {
  N <- length(y)
  fk <- fft(y) / N # normalize Data
  fk <- 2 * fk[1:((length(fk) / 2) + 1)] # DC comp + half of positives
  freq <- (0:(N - 1)) * sampFreq / N
  freq <- freq[(1:(length(fk)))]
  return(data.frame(fur = fk, freq = freq))
}
