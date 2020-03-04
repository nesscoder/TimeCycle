#' scale Time Series
#'
#' scales Time Series values to be between 0 and 1
#'
#' @param TimeSeries double
#'
#' @return double
#' @export
#'
scaleTimeSeries <- function(TimeSeries){
  minVal <- min(TimeSeries)

  if(sum(TimeSeries-mean(TimeSeries) == 0) == length(TimeSeries)){
    output <- TimeSeries
  } else if(minVal < 0){
    output <- TimeSeries + abs(minVal)
  } else{
    output <- TimeSeries - minVal
  }

  maxVal <- max(output)

  if(maxVal == 0){
    output <- rep(0, length(output))
  } else{
    output <- output/maxVal
  }

  return(output)
}
