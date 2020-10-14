#' Computes the Moving Average By Row in a \code{data.frame} of Time-Series
#'
#' Computes the moving average about each time-series in a \code{data.frame}.
#'
#' @param data a \code{data.frame} of \code{numeric} time-series expression values.
#'
#' @return a \code{data.frame} containing the smoothed \code{numeric} moving average time-series expression values by row.
#'
#' @seealso \code{\link{movingAverage}} for parameter definitions
#'
#' @export
#'
movingAverageDF <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)

  #smooth the data with an moving average filter
  dataDenoise <- apply(data, 1, function(TS) {
    movingAverage(x = as.numeric(TS), n = 3, centered = T)
  })

  output <- as.data.frame(t(as.data.frame(dataDenoise)))
  rownames(output) <- sampleNames
  colnames(output) <- colNames

  return(output)
}
