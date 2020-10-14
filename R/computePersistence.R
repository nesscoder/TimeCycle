#' Computes Persistence Scores For a Data.Frame of Time-Series Across Multiple Lags
#'
#' Takes a \code{data.frame} of numeric gene expression over time (genes X ZT times) and computes the persistence score using \code{\link{getPersistence}}.
#' For a given gene, each lag (min to max) is used to transform the expression into a 3-D embedded space via time-delay embedding.
#' A non-linear dimension reduction technique (laplacian eigenmaps) is used to transfrom the 3-D embedding to a 2-D embedding.
#' Finally, the persistence score of the 2-D embedding is calculated via persistence homology.
#' The median persistence score across all lags (min to max) for each gene is returned as a numeric vector.
#' For more details see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param minLag a \code{numeric} specifying the min lag to check in the 3-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 3-D embedding. Default is \code{5}.
#' @param cores a \code{numeric} specifying the number of parallel cores to use. Default number of cores is \code{parallel::detectedCores() - 2}.
#' @param laplacian a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}.
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
#'      \item \code{\link{getPersistence}} for use on a single gene expression time-series.
#'}
#'
#' @return a \code{vector} of the median persistence score across lags (minLag to maxLag) for each gene in data
#' @export
#'
computePersistence <- function(data, minLag = 2, maxLag = 5, cores = parallel::detectCores() - 2, laplacian = T){
  vect <- as.list(minLag:maxLag)
  #compute PS at each lag
  output <- parallel::mclapply(vect, mc.cores = cores, function(lag){
    perTSoutput <- apply(data,1,function(TS){
      return(getPersistence(t(as.matrix(TS)), lag = lag, laplacian = laplacian))
    })
    perTSoutput <- as.data.frame(perTSoutput)
  })

  #save persistence score at each lag
  PSatEachLag <- do.call(cbind,output)

  return(apply(PSatEachLag,1, stats::median))
}
