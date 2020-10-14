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
#'      \item Fasy, Brittany & Kim, Jisu & Lecci, Fabrizio & Maria, Clément. (2014). "Introduction to the R package TDA".
#'      \item Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }.
#'      \item Morozov D (2007). "Dionysus, a C++ library for computing persistent homology". \url{ http://www.mrzv.org/software/dionysus/ }
#'      \item Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
#'      \item Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.
#'    }
#' }
#' @seealso
#' \itemize{
#'      \item \code{\link[TDA]{ripsDiag}} for Persistence Homology calculation.
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
  persistence <- TDA::ripsDiag(X = embeddingLap, maxdimension = 1, maxscale = maxScale, library = "GUDHI", location = TRUE, printProgress = FALSE)

  # create a vectors loops formed by complex
  loops <- which(persistence$diagram[, 1] == 1)

  if (length(loops) > 0) {
    diff <- persistence$diagram[loops, 3] - persistence$diagram[loops, 2]

    # select persistence loops with max diff between Birth and Death
    maxdiff <- which.max(diff)
    birth <- persistence$diagram[loops, 2][maxdiff]
    death <- persistence$diagram[loops, 3][maxdiff]
    maxPersit <- unlist(unname(death - birth))
    return(maxPersit)
  } else {
    return(0)
  }
}
