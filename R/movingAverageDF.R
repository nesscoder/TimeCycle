#' movingAverageDF
#'
#' Computes the moving average across TimeSeries in a data frame
#'
#'
#' @param data double DataFrame /Rows = Genes, Column = Sample ZT/
#'
#' @return double DataFrame /Rows = Genes, Column = Sample ZT/
#' @export
#'
movingAverageDF <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)

  #smooth the data with an moving average filter
  dataDenoise <- apply(data, 1, function(TS) {
    movingAverage(x = as.numeric(TS), n = round(length(TS)/8), centered = T)
  })

  output <- as.data.frame(t(as.data.frame(dataDenoise)))
  rownames(output) <- sampleNames
  colnames(output) <- colNames

  return(output)
}
