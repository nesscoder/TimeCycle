#' getPersistence
#'
#' Returns the Max persistence score, returns 0 if no persistenec score exists
#'
#' @param TimeSeries double
#' @param lag int
#' @param laplacian boolean
#'
#' @return
#' @export
#'
getPersistence <- function(TimeSeries, lag, laplacian = T){

  #embed the Points
  embedding <- tryCatch({
    buildTakens_ndim(x = TimeSeries, dim = 3, delay = lag)
  }
  , warning = function(warning_condition) {
  }, error = function(error_condition) {
    #print("Insufficient observations for the requested embedding")
    return(NULL)
  }, finally={
  })

  if(is.null(embedding)){
    #print("Insufficient observations for the requested embedding, setting Persistence Score to 0")
    return(0)
  }

  #compute the Laplacian Embedding about the merged embedding
  embeddingLap <- tryCatch({
    embeddingLap <- computeLaplacianEmbedding(embedding)
  }
  , warning = function(warning_condition) {
    "warning"
  }, error = function(error_condition) {
    #print("Unable to Compute Eigen Vector, embedding colapsed")
    return(NULL)
  }, finally={
  })

  if(is.null(embeddingLap)){
    return(0) # setting Persistence Score to 0
  }

  #compute Rips Complex on the laplacian matrix
  #scale of Persistence to Check up to in persistence
  persistence <- TDA::ripsDiag(X = embeddingLap, maxdimension = 1, maxscale = 2, library = "GUDHI", location = TRUE, printProgress = FALSE)

  #create a vectors loops formed by complex
  loops <- which(persistence$diagram[,1] == 1)

  if(length(loops) > 0){
    diff <- persistence$diagram[loops,3] - persistence$diagram[loops,2]

    #select persistence loops with max diff between Birth and Death
    maxdiff <- which.max(diff)
    birth <- persistence$diagram[loops,2][maxdiff]
    death <- persistence$diagram[loops,3][maxdiff]
    maxPersit <- unlist(unname(death-birth))
    return(maxPersit)
  }else{
    return(0)
  }
}
