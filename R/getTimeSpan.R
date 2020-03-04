#' Get Differene Between Start Time and End Time
#'
#' Converts System Time Differences to
#' hh:mm:ss format
#'
#' @param start double
#' @param end double
#'
#' @return double
#' @export
#'
#' @examples
#' getTimeSpan(Sys.time(),Sys.time()+600)
getTimeSpan <- function(start, end) {
  dsec <- round(abs(as.numeric(difftime(end, start, units = "secs"))))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}
