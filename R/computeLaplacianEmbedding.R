#' Computes the Laplacian Eigenmap For Dimension Reduction
#'
#' Takes a  n-D Takens' Embedding and returns the 2-D Laplacian Eigenmap.
#'
#' @param takens a \code{data.frame} of the n-dimensional Takens' embedding. Columns defined from (1-D to n-D). TimeCycle defaults to \code{n = 3}.
#'
#' @return a \code{data.frame} of the 2-dimensional Takens' embedding. Columns defined by eigen vectors with last two non-trivial eigen values.
#'
#' @seealso
#' \itemize{
#'      \item \code{\link{getPersistence}} for usage
#'      \item \code{\link{matrixLaplacian}} for defining laplacian eigenmaps.
#'}
#'
#'
#' @export
#'
computeLaplacianEmbedding <- function(takens) {
  A <- as.matrix(stats::dist(x = takens, method = "euclidian"))
  A <- 1 - A / max(A) # normalize

  Z <- matrixLaplacian(A)
  # select last two non trivial eigen values
  eigen1 <- dim(Z$eigenvector)[2] - 2
  eigen2 <- dim(Z$eigenvector)[2] - 1
  return(Z$eigenvector[, c(eigen1, eigen2)])
}
