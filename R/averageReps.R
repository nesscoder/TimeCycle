#' averageReps
#'
#' averages expr values within replicate groups
#'
#' @param geneExpr double
#' @param Reps int
#'
#' @return double
#' @export
#'
#' @examples
#'
#' geneExpr <- c(1, 5, 3, 4, 5, 8)
#' reps <-  c(2, 3, 1)
#' avergeReps(geneExpr = geneExpr, Reps = reps)
avergeReps <- function(geneExpr,Reps){
  splitBy <- unlist(sapply(1:length(Reps), function(seq){
    rep(seq,Reps[seq])
  }))
  as.vector(unlist(lapply(split(geneExpr,splitBy),mean)))
}
