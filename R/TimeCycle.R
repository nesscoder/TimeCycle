#' Main TimeCycle Function
#'
#' Main function of the \pkg{TimeCycle} package used for detecting rhythmic signals in time-series gene expression sets.
#' For additional help with parameter selection, see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabel a \code{vector} defining the number of replicates at each time points.
#' @param resamplings a \code{numeric} specifying the number of resamplings to use in the null-distribution. Default is \code{10000}.
#' @param minLag a \code{numeric} specifying the min lag to check in the 3-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 3-D embedding. Default is  \code{5}.
#' @param cores a \code{numeric} specifying the number of parallel cores to use. Default number of cores is \code{parallel::detectedCores() - 2}.
#' @param period a \code{numeric} specifying the period of interest in hours for rhythm detection. Default is \code{24}.
#' @param laplacian  a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}
#' @param linearTrend  a \code{logical} scalar. Should TimeCycle Prioritize detecting linear trending signals? Default \code{FALSE}. Not recommended to change from default \code{FALSE} - will increases false positives rate. See vignette("TimeCycle") for more details.
#'
#' @references{
#' \itemize{
#'    \item{A pre-print version of TimeCycle is available on BioRxiv at \url{https://doi.org/10.1101/2020.11.19.389981}}
#'    }
#'    \subsection{TDA Package References}{
#'    \itemize{
#'    \item{Fasy, Brittany & Kim, Jisu & Lecci, Fabrizio & Maria, Clément. (2014). "Introduction to the R package TDA".}
#'    \item{Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }}
#'    \item{Morozov D (2007). "Dionysus, a C++ library for computing persistent homology". \url{ http://www.mrzv.org/software/dionysus/ }}
#'    \item{Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.}
#'    \item{Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.}
#'    }
#'
#'
#'    }
#'}
#'
#'
#'
#' @return a tidy \code{data.frame} of processed results for each gene:
#' \tabular{lcccccc}{
#'  \strong{sampleName} \tab  \strong{perScore} \tab \strong{pVals} \tab \strong{pVals.adj} \tab \strong{Period.in.Hours} \tab \strong{Amp} \tab \strong{Phase.In.Hours} \cr
#'  the gene name \tab  the median persistence score across all lags (min to max)  \tab raw empirical \emph{p}-value \tab FDR adjusted \emph{p}-value \tab period (h) \tab amplitude \tab phase (h) \cr
#' }
#'
#' @examples
#' # use built in zhang2014 data set sampled every
#' # 2 hours for 48 hours (i.e. 24 time points with 1 replicate each).
#' # Search for genes with period of 24 hours.
#'
#' #set seed for reproducibility with random variables in example usage
#' #> set.seed(1234)
#'
#' #> TimeCycleResults <- TimeCycle(data = zhang2014,
#' #>                               repLabel = rep(1,24),
#' #>                               period = 24)
#'
#' # Check number of genes with FDR < 0.05 and period between 22 to 26 hours.
#' #> library(tidyverse)
#'
#' #> TimeCycleResults %>%
#' #>    filter(22 < Period.in.Hours & Period.in.Hours < 26) %>%
#' #>    filter(pVals.adj < 0.05) %>%
#' #>    glimpse()
#'
#'@export
TimeCycle <- function(data,
                      repLabel,
                      resamplings = 10000,
                      minLag = 2,
                      maxLag = 5,
                      period = 24,
                      cores = parallel::detectCores()-2,
                      laplacian = T,
                      linearTrend = F
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
  dataAvg <- getRepAvgedDataFrame(data = data, repLabel = repLabel)

  #impute Missing Values
  dataAvg <- imputeMissingData(dataAvg)

  #center the data and merge data back together
  centeredData <- meanCenter(dataAvg)

  #detrend the data and ONLY USE FOR Compute the period
  dataDetrend <- detrend(centeredData)

  print("Computing Periods")
  #Compute Periods of Replicate 1 Time Series from the detrended data - does not break FFT assumptions
  periods <- periodFinder(movingAverageDF(dataDetrend))

  #scale Data Between 0 and 1 across all genes
  #allow us to compare all genes at once to the null distribution rather than just one at a time
  dataScaled  <- as.data.frame(t(apply(centeredData, 1, scaleTimeSeries)))
  colnames(dataScaled) <- colnames(centeredData)

  #Smooth data with Autocorrelation
  preProcessedData <- preprocess_acf(dataScaled, period, linearTrend = linearTrend)

  ##----------------------- pre-process NUll Distribution Data -----------------------
  print("Pre-Processing Null Distribution")

  #create the Null resampling of the Data
  resampledDataForNull <- nullResampling(dataScaled, numExperiments = resamplings)

  #Smooth -> Mean Center
  resampledprocessedData <- preprocess_acf(resampledDataForNull, period, linearTrend = linearTrend)

  ##------------------------ Compute the NUll Distribution  ------------------------
  print("Computing Null Distribution")

  nullDist <- computePersistence(resampledprocessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

  #Order the Null Distribtion to Speed up Ranking
  nullDistOrder <- as.numeric(nullDist[order(as.vector(nullDist))])

  ##--------------------Calculate the Peristence Score for Each Sample ------------------------
  print("Computing Persistence Scores")
  dataPS <- computePersistence(preProcessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

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
  pVals <- unlist(pVals)/as.numeric(resamplings)
  pVals.adj <- stats::p.adjust(pVals,method = 'fdr')
  names(pVals.adj) <- sampleNames

  endTime <- Sys.time()
  runTime <- getTimeSpan(startTime,endTime)

  results <- data.frame(sampleNames,perScore,pVals,pVals.adj,t(periods[,]))

  print(paste0("TimeCycle Completed"))
  print(paste0("Analysis Time: ", runTime))

  return(results)
}
