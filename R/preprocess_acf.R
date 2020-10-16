#' Computes Time-Series Autocovariance Function
#'
#' Computing the autocovariance function from a \code{data.frame} of time-series. The resulting autocovariance
#' function is smoothed using a moving average defined by the period and sampling scheme.
#' Returns the smoothed time-series as a \code{data.frame}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param period a \code{numeric} specifying the period of interest in hours for rhythm detection. Default is \code{24}.
#' @param linearTrend  a \code{logical} scalar. Should TimeCycle Prioritize detecting linear trending signals? Default \code{FALSE}. Not recommended to change from default \code{FALSE} - will increases false positives rate. See vignette("TimeCycle") for more details.
#'
#' @return a smoothed \code{data.frame} of \code{numeric} gene expression covariance over time (row = genes \emph{x} col = ZT times).
#'
#' @seealso \code{\link{meanCenter}}, \code{\link{scaleTimeSeries}}
#'
#' @export
#'
#'
preprocess_acf <- function(data, period = 24, linearTrend = F){
  data <- as.data.frame(data)
  xVals <- as.numeric(colnames(data))
  len <- max(xVals)
  interval <- mean(diff(xVals))
  maxAcfLag <- min(dim(data)[2]-1, 2*period/interval)
  output <- apply(data, 1, function(ts){
    ts <- unlist(unname(ts))
    corr <- stats::acf(as.vector(ts), lag = maxAcfLag, plot = F, type = "covariance")
    corr <- as.vector(corr$acf)
    corr <- scale(corr)

    if(linearTrend){
      # better for Linear Trends, increases Sigmoid false positive
      toCheck <- diff(as.vector(movingAverage(corr, n = period/ceiling(len/period)/interval+1, centered =  T)))
    } else{
      toCheck <- as.vector(movingAverage(corr, n = period/ceiling(len/period)/interval+1, centered =  T))
    }

  })

  output <- as.data.frame(t(output))
  #return the smoothed, mean centered data
  return(output)
}
