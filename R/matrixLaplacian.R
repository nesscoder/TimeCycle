#' Compute Matrix Laplacian About Takens Embedding
#'
#' @param A int
#' @param plot2D boolean
#' @param plot3D boolean
#'
#' @return
#' @export
#'
matrixLaplacian <- function(A, plot2D = FALSE, plot3D = FALSE) {
  B <- A
  D <- matrix(0, nrow = dim(A)[1], ncol = dim(A)[1])
  diag(D) <- B %*% rep(1, dim(A)[1])
  diag(D) <- 1 / sqrt(diag(D))
  Q <- D %*% B %*% D
  N <- diag(1, dim(Q)[1]) - Q
  Eigen <- eigen(N)$vectors
  if (plot2D) {
    plot(Eigen[, 1], Eigen[, 2], xlab = "", ylab = "")
  }
  if (plot3D) {
    scatterplot3d::scatterplot3d(Eigen[, 1], Eigen[, 2], Eigen[, 3],
      xlab = "",
      ylab = "", zlab = ""
    )
  }
  object <- list(LaplacianMatrix = N, eigenvector = Eigen)
  return(object)
}
