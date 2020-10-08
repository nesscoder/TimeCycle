#' smoothData.sgolay
#'
#' Smooth Data Using Savitzky–Golay filter
#'
#' @param data double DataFrame /Rows = Genes, Column = Sample ZT/
#'
#' @return double Savitzky-Golay Filtered by Row
#' @export
#'
smoothData.sgolay <- function(data) {
  sampleNames <- rownames(data)
  colNames <- colnames(data)

  # smooth the data with an Sgolay filter
  dataDenoise <- apply(data, 1, function(TS) {
    signal::sgolayfilt(x = as.numeric(TS))
  })

  output <- as.data.frame(t(as.data.frame(dataDenoise)))
  rownames(output) <- sampleNames
  colnames(output) <- colNames

  return(output)
}
