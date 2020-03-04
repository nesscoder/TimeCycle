#' Preprocess
#'
#' Preprocesses the data by smoothing with an sgolay filter, mean centering, and scaling
#'
#' @param data double dataframe /Rows = Genes, Column = Sample ZT/
#'
#' @return double preprocessed dataframe /Rows = Genes, Column = Sample ZT/
#' @export
#'
preprocess <- function(data){
  data <- as.data.frame(data)

  #Smooth Data with sgolay filter
  ###THINK ABOUT DOES IT MAKE SENSE TO SMOOTH OR NOT?
  ###MORE REPS SAMPLE POINTS YES, LESS REPS SAMPLE POINTS NO?
  dataSmooth <- smoothData.sgolay(data)

  #mean center the data
  dataCenter <- meanCenter(dataSmooth)

  #scale data between 0 and 1
  dataScaled  <- as.data.frame(t(apply(dataCenter, 1, scaleTimeSeries)))

  #return the detrended, smoothed, mean centered data
  return(dataScaled)
}
