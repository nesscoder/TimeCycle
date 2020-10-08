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
#' @param period int
#' @param laplacian boolean
#'
#' @return
#' @export
#'
TimeCycle <- function(data,
                      repLabel,
                      experiments = 10000,
                      minLag = 2,
                      maxLag = 5,
                      period = 24,
                      cores = parallel::detectCores()-2,
                      laplacian = T
){

  ## -----------------------------------pre-process Data ----------------------------

  cat("
      ########################################################################################
      ###      ████████ ██ ███    ███ ███████  ██████ ██    ██  ██████ ██      ███████     ###
      ###         ██    ██ ████  ████ ██      ██       ██  ██  ██      ██      ██          ###
      ###         ██    ██ ██ ████ ██ █████   ██        ████   ██      ██      █████       ###
      ###         ██    ██ ██  ██  ██ ██      ██         ██    ██      ██      ██          ###
      ###         ██    ██ ██      ██ ███████  ██████    ██     ██████ ███████ ███████     ###
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
  dataDetrend <- detrend(centeredData)

  print("Computing Periods")
  #Compute Periods of Replicate 1 Time Series from the detrended data - does not break FFT assumptions
  periods <- periodFinder(movingAverageDF(dataDetrend),plotFit = F)

  #scale Data Between 0 and 1 across all genes
  #allow us to compare all genes at once to the null distribution rather than just one at a time
  dataScaled  <- as.data.frame(t(apply(centeredData, 1, scaleTimeSeries)))
  colnames(dataScaled) <- colnames(centeredData)
  #Smooth -> Mean Center
  dataScaled  <- dataScaled #store as global variable to be accessed by the user
  preProcessedData <- preprocess_acf(dataScaled, period)

  ##----------------------- pre-process NUll Distribution Data -----------------------
  print("Pre-Processing Null Distribution")

  #create the Null resampling of the Data
  resampledDataForNull <- nullResampling(dataScaled, numExperiments = experiments)

  #Smooth -> Mean Center
  resampledprocessedData <- preprocess_acf(resampledDataForNull, period)

  ##------------------------ Compute the NUll Distribution  ------------------------
  print("Computing Null Distribution")

  nullDist <- computePersistence(resampledprocessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

  #Order the Null Distribtion to Speed up Ranking
  nullDistOrder <- as.numeric(nullDist[order(as.vector(nullDist))])

  ##--------------------Calculate the Peristence Score for Each Sample ------------------------
  print("Computing Persistence Scores")
  dataPS <<- computePersistence(preProcessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

  ##-----------------------    Get pVals for Data ------------------------
  print("Calculating p-values")
  perScore <- as.numeric(dataPS)
  sampleNames <- rownames(preProcessedData)
  vect <- 1:length(perScore)

  #Compute Emperical pVal
  pVals <- parallel::mclapply(X = vect, mc.cores = cores, function(x){
    return(rank(-c(perScore[x],nullDistOrder),ties.method = "max")[1])
  })

  #get FDR Adjusted pVal
  pVals <- unlist(pVals)/as.numeric(experiments)
  pVals.adj <- stats::p.adjust(pVals,method = 'fdr')
  names(pVals.adj) <- sampleNames

  endTime <- Sys.time()
  runTime <- getTimeSpan(startTime,endTime)

  results <- data.frame(sampleNames,perScore,pVals,pVals.adj,t(periods[,]))

  print(paste0("TimeCycle Completed"))
  print(paste0("Analysis Time: ", runTime))

  return(results)
}
