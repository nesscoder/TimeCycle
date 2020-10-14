#' Mean Centers a \code{data.frame}
#'
#' Mean centers a \code{data.frame} by row.
#'
#' @param  df a \code{data.frame} of \code{numerics}.
#'
#' @return a mean center \code{data.frame} of \code{numerics} by row.
#' @export
#'
meanCenter <- function(df){
  as.data.frame(t(apply(df, 1, function(y) y-mean(y))))
}
