#' Detrend
#'
#' Detrend Data by fiting a linear model and removing the linear trend
#'
#' @param data double dataframe /Rows = Genes, Column = Sample ZT/
#'
#' @return double detrended dataframe /Rows = Genes, Column = Sample ZT/
#' @export
#'
detrend <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)
  #get numeric timepoints values from column Names
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #Use timeseries to compute the trend fit
  dataDetrend <- apply(data, 1, function(TS) {
    fit <- lm(TS~xVals)
    y <- zapsmall(TS - (xVals*fit$coefficients[2] + fit$coefficients[1]),10)
  })

  #return the detrended data as a data.frame
  dataDetrend <- as.data.frame(t(as.data.frame(dataDetrend)))
  rownames(dataDetrend) <- sampleNames
  colnames(dataDetrend) <- colNames

  return(dataDetrend)
}
