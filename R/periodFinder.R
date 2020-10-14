#' Estimates the period, amp, and phase of a \code{data.frame} of time-series
#'
#' Estimates the period, amp, and phase of a \code{data.frame} of time-series via fast Fourier transform (FFT).
#' A model fit for the time-series is generated using the first 3 harmonics. The period, amp, and phase are
#' computed based on the aggregate fit.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @return  a \code{data.frame} of \code{numeric} period, amplitude, and phase estimates for each gene.
#' @seealso
#' \itemize{
#'      \item \code{\link{getFFT}} for FFT calculation.
#'      \item \code{\link{detrend}} for linear detrending of time-series prior to periodFinder processing.
#'      \item \code{\link{movingAverageDF}} for smoothing of time-series prior to periodFinder processing.
#'}
#'
#'
#' @export
#'
periodFinder <- function(data){

  colNames <- colnames(data)

  #get Time Points from Colnames
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #get difference between Time Points
  avgDiff <- mean(diff(xVals))

  dataPeriod <- apply(data, 1, function(TS) {

    #compute the FFT on the Time series
    fft <- getFFT(TS,1/avgDiff)

    #create variable to hold FFT fits with first 3 harmonics
    fit <- list()
    for(i in 1:4){

      #Get period
      per <- (1/fft$freq[i])

      #Get Amplitude
      amp <- Mod(fft$fur)[i]

      #Get Phase
      phase <- Arg(fft$fur)[i]

      #get Cos Wave
      tp <- seq(0,diff(range(xVals)),0.1)
      harmonicFit <- amp*cos(2*pi*(tp/per)+phase)

      #save the harmonic cos fit
      fit[[i]] <- as.vector(harmonicFit)
    }

    #compute Signal using first 3 FFT harmonics
    fftFullfit <- rowSums(do.call(cbind,fit))

    #get local Max and local Min of signal
    localMax <- (which(diff(sign(diff(fftFullfit)))==-2)+1)
    localMin <- (which(diff(sign(diff(fftFullfit)))== 2)+1)
    minMax <- c(localMax,localMin)
    minMaxOrdered <- minMax[order(minMax)]
    deriv <- diff(fftFullfit)

    #search for outlier points identified as min and maxes
    pos <- c(1,minMaxOrdered)
    pat <- rep(seq_along(pos), times=diff(c(pos, length(fftFullfit))))
    z <- split(abs(deriv)/max(abs(deriv)), pat)

    #for each point compute the area under the derivative curve
    areas <- as.vector(unlist(lapply(z,sum)))/sum(abs(deriv)/max(abs(deriv)))

    #if normalized area is less than 0.05 of total, consider outlier
    outliers <- which(areas[2:length(areas)] < 0.05)
    #if normalized area is more than 0.05 of total, keep point
    keep <- which(areas[2:length(areas)] > 0.05)

    for(i in 1:length(keep)){
      if(length(outliers) > 0){
        toAvg <- which(outliers < keep[i]) #which outlier points are less than points to keep
        if(length(toAvg) > 0 ){ #check if there are points to average
          start <- outliers[toAvg][1] #start from outlier point less than keep point
          stop <- keep[i] #upto the first keep point
          selPoints <- start:stop
          minMaxOrdered[keep[i]] <- mean(minMaxOrdered[selPoints]) #average points together
          outliers <- outliers[-toAvg] #remove the points from outlier
        }
      }

    }
    minMaxOrdered <- minMaxOrdered[keep]

    #Get overall period
    if(length(localMax)+length(localMin) < 2){
      #if only one max and min, the period is the sampling length
      per <- abs(diff(range(xVals)))
    } else{
      #else the period is the avg diff between all local maxima and minima
      per <- round(mean(2*abs(diff(minMaxOrdered)))*0.1,2)
    }

    #Get overall Amplitude
    #mean of the absolute value of the mean center signal. taking all local min and max into account
    amp <- mean(abs(fftFullfit[minMaxOrdered] - mean(fftFullfit)))

    #Get overall Phase
    #time to first localMax. ie. how far is the cos wave shifted from zero
    if(length(localMax) < 1){
      phase <-  0
    }else{
      phase <- localMax[1]
    }
    phaseInHours <- (phase*0.1)%%per

    return(round(c(per, amp, phaseInHours),2))
  })

  dataPeriod <- as.data.frame(dataPeriod)
  rownames(dataPeriod) <- c("Period.in.Hours", "Amp", "Phase.in.Hours")
  return(dataPeriod)
}
