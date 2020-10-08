#' computeLaplacianEmbedding
#'
#' Compute the Laplacian Embedding about the merged Takens embedding
#'
#' @param takens double
#'
#' @return double
#' @export
#'
computeLaplacianEmbedding <- function(takens) {
  A <- as.matrix(dist(x = takens, method = "euclidian"))
  A <- 1 - A / max(A) # comment out if unnorms

  Z <- matrixLaplacian(A, plot2D = FALSE, plot3D = FALSE)
  # select last two non trivial eigen values
  eigen1 <- dim(Z$eigenvector)[2] - 2
  eigen2 <- dim(Z$eigenvector)[2] - 1
  return(Z$eigenvector[, c(eigen1, eigen2)])
}
