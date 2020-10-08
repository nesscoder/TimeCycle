#' Null Resampling
#'
#' Generates a null resampling as defined by rhe difference between consecutive points across all time series
#'
#' @param data double dataframe /Rows = Genes, Column = Sample ZT/
#' @param numExperiments int number of Resamplings to perform
#'
#' @return double dataframe /Rows = Genes, Column = Sample ZT/
#' @export
#'
nullResampling <- function(data, numExperiments){

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
