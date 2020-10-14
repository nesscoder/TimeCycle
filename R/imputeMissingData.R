#' Imputates Missing Time-Points in Data
#'
#' Imputes \code{numeric} values for time-points with an \code{NA} by computing the linear path between missing points
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times) with missing values.
#'
#' @return a imputed \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @export
#'
#' @seealso \code{\link[imputeTS]{na_interpolation}} for imputation procedure.
#'
#' @examples
#' a <- c(10, 12, 14, NA, 18)
#' b <- c(1, 2, NA, NA, 5)
#' data <- t(data.frame(a, b))
#' imputeMissingData(data)
imputeMissingData <- function(data) {
  colNames <- colnames(data)
  rowNames <- rownames(data)

  out <- t(apply(data, 1, function(geneExpr) {
    x <- as.numeric(geneExpr)
    imputed <- imputeTS::na_interpolation(x)
  }))

  colnames(out) <- colNames
  rownames(out) <- rowNames

  return(out)
}
