#' computePersistence2.0
#'
#' Returns the Lag with the Top Persistence Score
#'
#' @param data int
#' @param minLag int
#' @param maxLag int
#' @param cores int
#' @param laplacian boolean
#'
#' @return int
#' @export
#'
computePersistence <- function(data, minLag = 2, maxLag = 5, cores = 1, laplacian = T){
  vect <- as.list(minLag:maxLag)
  .xVals <<- as.numeric(unlist(regmatches(colnames(data), gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colnames(data)))))
  #compute PS at each lag
  output <- mclapply(vect, mc.cores = cores, function(lag){
    perTSoutput <- apply(data,1,function(TS){
      return(getPersistence(t(as.matrix(TS)), lag = lag, laplacian = laplacian))
    })
    perTSoutput <- as.data.frame(perTSoutput)
  })

  #save persistence score at each lag
  PSatEachLag <- do.call(cbind,output)

  return(apply(PSatEachLag,1, median))
}
