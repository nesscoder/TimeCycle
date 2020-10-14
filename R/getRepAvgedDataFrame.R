#' Handles Sample Replicates in a \code{data.frame} of Time-Series
#'
#' Averages expression values for a single time-series across replicate groups in a \code{data.frame}.
#'
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabel a \code{vector} defining the number of replicates at each time point.
#'
#' @return a \code{data.frame} of average  \code{numeric} time-series expression values by replicate time-points for each gene.
#'
#' @seealso \code{\link{averageReps}}
#'
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
