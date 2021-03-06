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
#'    \item{The TimeCycle Manuscript is available from Bioinformatics at \url{https://doi.org/10.1093/bioinformatics/btab476}}
#'    }
#'    \subsection{TDA Package References}{
#'    \itemize{
#'      \item Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. (2018). "TDAstats: R pipeline for computing persistent homology in topological data analysis." \emph{Journal of Open Source Software}. 2018; 3(28): 860. doi:\href{https://doi.org/10.21105/joss.00860}{[10.21105/joss.00860]}
#'      \item Bauer U. (2019). "Ripser: Efficient computation of Vietoris-Rips persistence barcodes." \emph{arXiv}: 1908.02518.
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
