#' Computes Persistence Scores For a Single Time-Series Across a Single Lag
#'
#' Takes a \code{vector} of numeric gene expression over time and computes the persistence score.
#' The specified lag is used to transform the expression into a 3-D embedded space via time-delay embedding.
#' A non-linear dimension reduction technique (laplacian eigenmaps) is used to transfrom the 3-D embedding to a 2-D embedding.
#' Finally, the persistence score of the 2-D embedding is calculated via persistence homology.
#' Returns the Max persistence score, returns 0 if no persistence score exists.
#' For more details see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param timeSeries a \code{vector} of \code{numeric} time-series expression values.
#' @param lag a \code{numeric} specifying the Lag to use for in the 3-D time delayed embedding.
#' @param laplacian a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}.
#'
#' @return the max persistence score at the specified lag, returns 0 if no persistence score exists.
#'
#' @references{
#'    \itemize{
#'      \item Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. (2018). "TDAstats: R pipeline for computing persistent homology in topological data analysis." \emph{Journal of Open Source Software}. 2018; 3(28): 860. doi:\href{https://doi.org/10.21105/joss.00860}{[10.21105/joss.00860]}
#'      \item Bauer U. (2019). "Ripser: Efficient computation of Vietoris-Rips persistence barcodes." \emph{arXiv}: 1908.02518.
#'    }
#' }
#' @seealso
#' \itemize{
#'      \item \code{\link[TDAstats]{calculate_homology}} for Persistence Homology calculation.
#'      \item \code{\link{buildTakens_ndim}} for for generating time-delay embedding.
#'      \item \code{\link{computeLaplacianEmbedding}} for 3-D to 2-D laplacian eigenmaps dimension reduction.
#'      \item \code{\link{computePersistence}} for use parallelized function for a \code{data.frame} of gene expression.
#'    }
#' @export
#'
getPersistence <- function(timeSeries, lag, laplacian = T) {


  # embed the Points
  embedding <- tryCatch(
    {
      if (laplacian) {
        buildTakens_ndim(x = as.vector(unlist(timeSeries)), dim = 3, delay = lag)
      } else {
        buildTakens_ndim(x = as.vector(unlist(timeSeries)), dim = 2, delay = lag)
      }
    },
    warning = function(warning_condition) {
    },
    error = function(error_condition) {
      # print("Insufficient observations for the requested embedding")
      return(NULL)
    },
    finally = {
    }
  )


  if (is.null(embedding)) {
    # print("Insufficient observations for the requested embedding, setting Persistence Score to 0")
    return(0)
  }

  # compute the Laplacian Embedding about the merged embedding
  embeddingLap <- tryCatch(
    {
      if (laplacian) {
        embeddingLap <- computeLaplacianEmbedding(embedding)
      } else {
        embeddingLap <- embedding
      }
    },
    warning = function(warning_condition) {
      "warning"
    },
    error = function(error_condition) {
      # print("Unable to Compute Eigen Vector, embedding colapsed")
      return(NULL)
    },
    finally = {
    }
  )

  if (is.null(embeddingLap)) {
    return(0) # setting Persistence Score to 0
  }

  # compute Rips Complex on the laplacian matrix
  # scale of Persistence to Check up to in persistence
  maxScale <- ceiling(2 * max(abs(as.vector(unlist(timeSeries)))))
  maxScale <- max(maxScale, 1, na.rm = T)
  persistence <- TDAstats::calculate_homology(embeddingLap, dim = 1, threshold = maxScale)
  # create a vectors loops formed by complex
  loops <- which(persistence[, 1] == 1)

  if (length(loops) > 0) {

    # select persistence loops with max diff between Birth and Death
    diff <- persistence[loops, 3] - persistence[loops, 2]
    maxdiff <- which.max(diff)
    maxPersit <- diff[maxdiff]
    return(maxPersit)
  } else {
    return(0)
  }
}
