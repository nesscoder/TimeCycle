#' periodFinder
#'
#' Computes the period, amp, and phase of a dataframe of time series
#'
#' @param data double dataframe (row = genes, col = timepoints)
#' @param plotFit boolean
#'
#' @return
#' @export
#'
periodFinder <- function(data, plotFit = F){

  colNames <- colnames(data)

  #get Time Points from Colnames
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #get difference between Time Points
  avgDiff <- mean(diff(xVals))

  dataPeriod <- apply(data, 1, function(TS) {
    #Generate Periodogram of the signal

    periodogram <- TSA::periodogram(as.numeric(TS),log='no',plot=F)

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
    if(plotFit){
      par(mfrow = c(2,2))
      plot(fftFullfit)
      points(minMaxOrdered,fftFullfit[minMaxOrdered], col = "red", pch = 16)
      plot(deriv)
      abline(h = 0, col = "red")
      call <- c(abs(deriv)/max(abs(deriv)) < 0.05)
      plot(abs(deriv)/max(abs(deriv)), col = call+1, pch = 16)
      abline(h = 0.05, col = "red")
    }

    #search for outlier points identified as min and maxes
    pos <- c(1,minMaxOrdered)
    pat <- rep(seq_along(pos), times=diff(c(pos, length(fftFullfit))))
    z <- split(abs(deriv)/max(abs(deriv)), pat)

    #for each point compute the area under the derivative curve
    areas <- as.vector(unlist(lapply(z,sum)))/sum(abs(deriv)/max(abs(deriv)))

    if(plotFit){
      plot(areas[2:length(areas)], col = (areas[2:length(areas)] < 0.05)+1)
    }

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

    #plot fit to Signal
    tp <- seq(0,diff(range(xVals)),0.1)
    if(plotFit){
      par(mfrow = c(2,2))
      plot(xVals-min(xVals), TS,"o", xlim = c(0,diff(range(xVals))), pch = 16, main = paste0("Period: ", round(per,2)," | Amp: ", round(amp,2)," | Phase: ", round(phaseInHours,2)))
      abline(h = 0, col = "grey", lty = 2)
      for(i in seq(0,max(xVals),per)){
        abline(v = i, col = "grey", lty = 2)
      }

      #plot individual harmonics
      points(tp, fftFullfit, lwd=3, type = "l", lty=1, col="blue")
      for( i in seq_along(fit)){
        points(tp, fit[[i]], lwd=2, type = "l", lty=2, col=i)
      }

      #plot phase and amplitude
      if(phase != 0){
        shiftedMax <- tp[phase]
      }else{
        shiftedMax <- phase
      }
      segments(x0 = 0, y0 = amp, x1 = phaseInHours, y1 = amp, col = "orange", lwd = 3)
      segments(x0 = shiftedMax, y0 = 0, x1 = shiftedMax, y1 = amp, col = "green", lwd = 3)
      axis(side = 1,at = seq(0,max(xVals),5), labels = seq(0,max(xVals),5))

      #Plot Periodgram information
      plot(periodogram$freq,periodogram$spec, "o", main = "Periodigram", xlab = "freq", ylab = "Spectrum Intensity")
      plot(fft$freq,Mod(fft$fur), "o", main = "Amp Spec", xlab = "freq", ylab = "Amp")
      plot(fft$freq,Arg(fft$fur), "o", main = "Phase Spec", xlab = "freq", ylab = "Phase" )
    }

    return(round(c(per, amp, phaseInHours),2))
  })

  dataPeriod <- as.data.frame(dataPeriod)
  rownames(dataPeriod) <- c("Period.in.Hours", "Amp", "Phase.in.Hours")
  return(dataPeriod)
}
