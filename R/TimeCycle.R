#' TimeCycle
#'
#' Main Function
#'
#' @param data double
#' @param repLabel int
#' @param experiments int
#' @param minLag int
#' @param maxLag int
#' @param cores int
#'
#' @return
#' @export
#'
#' @examples
TimeCycle <- function(data, repLabel, experiments = 10000, minLag = 1, maxLag = 6, cores = detectCores()-2){

  ## -----------------------------------pre-process Data ----------------------------

  cat("
      ########################################################################################
      ### ████████╗ ██╗ ███╗   ███╗ ███████╗  ██████╗ ██╗   ██╗  ██████╗ ██╗      ███████╗ ###
      ### ╚══██╔══╝ ██║ ████╗ ████║ ██╔════╝ ██╔════╝ ╚██╗ ██╔╝ ██╔════╝ ██║      ██╔════╝ ###
      ###    ██║    ██║ ██╔████╔██║ █████╗   ██║       ╚████╔╝  ██║      ██║      █████╗   ###
      ###    ██║    ██║ ██║╚██╔╝██║ ██╔══╝   ██║        ╚██╔╝   ██║      ██║      ██╔══╝   ###
      ###    ██║    ██║ ██║ ╚═╝ ██║ ███████╗ ╚██████╗    ██║    ╚██████╗ ███████╗ ███████╗ ###
      ###    ╚═╝    ╚═╝ ╚═╝     ╚═╝ ╚══════╝  ╚═════╝    ╚═╝     ╚═════╝ ╚══════╝ ╚══════╝ ###
      ########################################################################################\n")

  #remove before launching
  # set.seed(123)

  startTime <- Sys.time()
  print("Starting TimeCycle")
  print("Pre-Processing Data")

  #reorder the Data By Replicate
  dataAvg <<- getRepAvgedDataFrame(data = data, repLabel = repLabel)

  #impute Missing Values
  dataAvg <- imputeMissingData(dataAvg)

  #center the data and merge data back together
  centeredData <- meanCenter(dataAvg)

  #detrend the data and ONLY USE FOR Compute the period
  dataDetrend <<- detrend(centeredData)

  #Compute Periods of Replicate 1 Time Series from the detrended data - does not break FFT assumptions
  periods <- periodFinder(movingAverageDF(dataDetrend),plotFit = F)

  #scale Data Between 0 and 1 across all genes
  #allow us to compare all genes at once to the null distribution rather than just one at a time
  dataScaled  <- as.data.frame(t(apply(centeredData, 1, scaleTimeSeries)))
  colnames(dataScaled) <- colnames(centeredData)
  #Smooth -> Mean Center
  dataScaled  <<- dataScaled #store as global variable to be accessed by the user
  preProcessedData <<- preprocess(dataScaled)

  ##----------------------- pre-process NUll Distribution Data -----------------------
  print("Pre-Processing Null Distribution")

  #create the Null resampling of the Data
  resampledDataForNull <<- nullResampling(dataScaled, numExperiments = experiments)

  #Smooth -> Mean Center
  resampledprocessedData <<- preprocess(resampledDataForNull)

  ##------------------------ Compute the NUll Distribution  ------------------------
  print("Computing Null Distribution")

  #Compute Null Distributions at Each Lag
  vect <- as.list(c(minLag:maxLag))

  nullDist <<- parallel::mclapply(vect, mc.cores = cores, function(lag){
    apply(resampledprocessedData,1,function(TS){
      getPersistence(t(as.matrix(TS)),as.numeric(lag),laplacian = T)
    })
  })

  #Order the Null Distribtion to Speed up Ranking
  nullDistOrder <-lapply(nullDist, function(nulld){
    nulld[order(as.vector(nulld))]
  })

  # PLOT the Null Distributions
  # This is a check and Not Needed for final Version
  #par(mfrow = c(1,3))
  #for(i in minLag:maxLag){
  #hist(nullDist[[i]], main = paste0("Hist of Null Dist Lag", i) , xlab = paste("Null Dist Lag", i), xlim = c(0,0.4), ylim = c(0, 5000))
  #}


  ##--------------------Calculate the Peristence Score for Each Sample ------------------------
  print("Computing Persistence Scores")
  dataPS <- computePersistence2.0(preProcessedData, minLag = minLag, maxLag = maxLag, cores = cores)

  ##-----------------------    Get pVals for Data ------------------------
  print("Calculating p-values")
  perScore <- as.numeric(dataPS[1,])
  perlag <- as.numeric(dataPS[2,])
  sampleNames <- rownames(preProcessedData)
  vect <- 1:length(perScore)

  #Compute Emperical pVal
  pVals <- parallel::mclapply(X = vect, mc.cores = cores, function(x){
    return(rank(-c(perScore[x],nullDistOrder[[perlag[x]]]))[1])
  })

  #get FDR Adjusted pVal
  pVals <- unlist(pVals)/as.numeric(experiments)
  pVals.adj <- p.adjust(pVals,method = 'fdr')
  names(pVals.adj) <- sampleNames

  endTime <- Sys.time()
  runTime <- getTimeSpan(startTime,endTime)

  results <- data.frame(sampleNames,perScore,perlag,pVals,pVals.adj,t(periods[,]))

  print(paste0("TimeCycle Completed"))
  print(paste0("Analysis Time: ", runTime))

  return(results)
}
