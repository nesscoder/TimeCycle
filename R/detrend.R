#' Linearly Detrends a \code{data.frame} of Gene Expression Time-Series
#'
#' Detrends a \code{data.frame} of gene expression over time (row = genes \emph{x} col = ZT times)
#' by fiting a linear model to each gene and removing the linear trend.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @return a detrended \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @seealso \code{\link{periodFinder}}
#'
#' @export
#'
detrend <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)
  #get numeric timepoints values from column Names
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #Use timeseries to compute the trend fit
  dataDetrend <- apply(data, 1, function(TS) {
    fit <- stats::lm(TS~xVals)
    y <- zapsmall(TS - (xVals*fit$coefficients[2] + fit$coefficients[1]),10)
  })

  #return the detrended data as a data.frame
  dataDetrend <- as.data.frame(t(as.data.frame(dataDetrend)))
  rownames(dataDetrend) <- sampleNames
  colnames(dataDetrend) <- colNames

  return(dataDetrend)
}
