#' computePersistence2.0
#'
#' Returns the Lag with the Top Persistence Score
#'
#' @param data int
#' @param minLag int
#' @param maxLag int
#' @param cores int
#'
#' @return int
#' @export
#'
computePersistence2.0 <- function(data, minLag = 1, maxLag = 6, cores = 1){
  vect <- as.list(minLag:maxLag)

  #compute PS at each lag
  output <- parallel::mclapply(vect, mc.cores = cores, function(lag){
    perTSoutput <- apply(data,1,function(TS){
      return(getPersistence(t(as.matrix(TS)), lag = lag, laplacian = T))
    })
    perTSoutput <- as.data.frame(perTSoutput)
  })

  #Select Lag with Top PS
  mergedOutput <- do.call(cbind,output)
  output <- apply(mergedOutput,1, function(PS){
    c(max(PS),which.max(PS))})

  return(unlist(output))
}
