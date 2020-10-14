#' Computes Laplacian Eigenmaping
#'
#' Computes the Laplacian matrix and eigenvectors of the n-D Takens' embedding.
#'
#' @param A a distance \code{matrix} of the n-dimensional Takens' embedding.
#'
#' @return a \code{list} of the laplacian matrix and eigenvectors.
#'
#' @seealso \code{\link{computeLaplacianEmbedding}}
#'
#' @export
#'
matrixLaplacian <- function(A) {
  B <- A
  D <- matrix(0, nrow = dim(A)[1], ncol = dim(A)[1])
  diag(D) <- B %*% rep(1, dim(A)[1])
  diag(D) <- 1 / sqrt(diag(D))
  Q <- D %*% B %*% D
  N <- diag(1, dim(Q)[1]) - Q
  Eigen <- eigen(N)$vectors

  object <- list(LaplacianMatrix = N, eigenvector = Eigen)
  return(object)
}
