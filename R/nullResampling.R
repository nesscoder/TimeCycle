#' Generates a \code{data.frame} of Resampled Time-Series for Computing the Null Distribution
#'
#' Generates a \code{data.frame} null resampling time-series as defined by the difference between consecutive points across all time-series.
#' For additional information about the null distribution used by TimeCycle, see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param numExperiments a \code{numeric} specifying the number of resampling to use in the null-distribution. Default is \code{10000}.
#'
#' @return a \code{data.frame} of \code{numeric} resampled time-series (row = numExperiments \emph{x} col = ZT times).
#'
#' @seealso \code{\link{resampleTimeSeries}}
#'
#' @export
#'
nullResampling <- function(data, numExperiments = 10000){

  #Get Average Value of Data across Replicate
  dataForResampling <- data
  colNames <- colnames(data)

  #get Resampled Data for Use in generating Null Distribution
  diffBetweenTP <- t(apply(dataForResampling, 1, function(geneName){
    #Get Difference Between Time Points
    diff(geneName)
  }))

  #Resample Differences
  resampledDataNull <- resampleTimeSeries(diffBetweenTP, numTP = dim(data)[2], numRsmps = numExperiments)
  colnames(resampledDataNull) <- colNames

  return(resampledDataNull)
}
