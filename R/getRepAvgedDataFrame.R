#' getRepAvgedDataFrame
#'
#' averages expr values within replicate groups by row in a dataframe
#'
#' @param data double dataframe /Rows = Genes, Column = Sample ZT/
#' @param repLabel int
#'
#' @return
#' @export
#'
getRepAvgedDataFrame <- function(data, repLabel) {

  # get Column Names
  colnames <- colnames(data)
  # remove replicate label from colnames if they exist
  colnames <- gsub(pattern = "_rep.", replacement = "", colnames)
  # get unique ZT time for each point
  colnames <- unique(as.numeric(unlist(regmatches(colnames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames)))))

  # get average of Replicate Labels
  output <- t(apply(data, MARGIN = 1, FUN = function(geneExprRow) {
    averageReps(geneExprRow, repLabel)
  }))
  output <- as.data.frame(output)
  colnames(output) <- colnames

  return(output)
}
