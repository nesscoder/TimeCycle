#' mean center the data by Row
#'
#' @param df double
#'
#' @return double
#' @export
#'
meanCenter <- function(df){
  as.data.frame(t(apply(df, 1, function(y) y-mean(y))))
}
