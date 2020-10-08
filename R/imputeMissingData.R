#' ImputeMissingData
#'
#' Imputes NA in the data frame by computing the linear path between missing points
#'
#' @param data df (Col = Genes, Row = TimePoints)
#'
#' @return
#' @export
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
