#' Preprocessing AutoCorrelation Function
#'
#' Take a dataframe of time-series as input and smooth oscillations by computing the autocorrelation function
#' and taking the moving average across time points according the defined period and sampling scheme.
#' Returns the smoothed time-series as a dataframe.
#'
#' @param data dataframe of time series
#' @param period the period of interest in your data (i.e. default = 24)
#'
#' @return
#' @export
#'
#'
preprocess_acf <- function(data, period = 24){
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

    #toCheck <- diff(as.vector(movingAverage(corr, n = length(ts)/4, centered =  T))) # better for Linear Trends, increases Sigmoid false positive
    toCheck <- as.vector(movingAverage(corr, n = period/ceiling(len/period)/interval+1, centered =  T))
  })

  output <- as.data.frame(t(output))
  #return the smoothed, mean centered data
  return(output)
}
