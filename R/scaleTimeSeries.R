#' Scales a Time-Series
#'
#' Scales a time-series between 0 and 1.
#'
#' @param TimeSeries a \code{vector} of \code{numeric} time-series expression values.
#'
#' @return a scaled \code{vector} of \code{numeric} time-series expression values between 0  and 1.
#'
#' @seealso \code{\link{preprocess_acf}}, \code{\link{meanCenter}}
#'
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
