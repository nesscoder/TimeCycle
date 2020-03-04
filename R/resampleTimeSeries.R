#' resampleTimeSeries
#'
#' Takes a dataframe of expr values and generates a resampling of the data used for parametric testing
#'
#' @param distTP double
#' @param numTP int
#' @param numRsmps int
#'
#' @return int dataframe
#' @export
#'
resampleTimeSeries <- function(distTP, numTP, numRsmps){
  dist <- as.numeric(as.vector(unlist(distTP)))
  return(as.data.frame(matrix(data = sample(dist, size = numTP*numRsmps, replace = T), ncol = numTP, nrow = numRsmps)))
}
