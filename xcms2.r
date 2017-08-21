#!/usr/bin/env Rscript
# @aly

# Packages required for all functions in this script to run
# Remember to install and load them in scripts that source this script
# suppressPackageStartupMessages({
#     library(e1071)
#     library(xcms)
#     library(nortest)
#     })

#Suppress scientific notation
options(scipen=999)

########################################################################################
# Parameters for centwave/getPeaks.
# Put peakwidths (PEAKWIDTHL) in ascending order. For quantitation, the order where the 
#  results are returned is in decending order (reverse).
# PEAKWIDTHG is used for getPeaks/pseudoPeakCalc baseline adjusted peak area.
########################################################################################
RTWINDOW <<- 30
#ppm value used if no ppm is specified, or if we do not provide an ROI to centWave
CWPPM <<- 25
PREFILTER <<- c(1, 1000)
SNTHRESH <<- 10
INTEGRATE <<- 2
MZDIFF <<- -10000.0
PEAKWIDTHL <<- list(c(2, 80), c(5, 80), c(10, 80), c(15, 80), c(20, 80))

PEAKWIDTHG <<- c(10, 80)
# GSTEP <<- 0.001

DEV <<- CWPPM*1e-6
MINCENTROIDS <<- 3

BASELINECOUNT <<- 3


########################################################################################
# Input: Mass, ppm
# Output: Mass range
# 
# This method calculates a mass range given a mass and a mass window (ppm).
# Numbers used are based on visual basic AutoFill
########################################################################################
pmMz <- function (x, ppm=10) {
    c(x - x * ppm * 0.000001, x + x * ppm * 0.000001)
}

########################################################################################
# Input: Retention time, retention time window (one side)
# Output: Retention time range
# 
# This method calculates a retention time window given a retention time and a retention
#  time window.
########################################################################################
pmRt <- function (x, rtwidth=10) {
    c(x - rtwidth, x + rtwidth)
}

########################################################################################
# Input: xcmsRaw object, retention time range
# Output: Scan number range
# 
# This method takes a retention time range, and converts it to a scan number range, by 
#  finding the closest retention time of the scan that is less than the minimum, and the
#  the closest retention time of the scan that is greater than the maximum, of the 
#  provided retention time range.
# This method is based on xcms source code. 
########################################################################################
rtRangeToScanRange <- function (xcmsRaw, rtrange) {
    rtrange <- range(rtrange)
    scanidx <- (xcmsRaw@scantime >= rtrange[1]) & (xcmsRaw@scantime <= rtrange[2])
    c(match(TRUE, scanidx), length(scanidx) - match(TRUE, rev(scanidx)))
}

########################################################################################
# Input: xcmsRaw object, acquisition number
# Output: Retention time closest to acquisition number
# 
# This method finds the retention time closest to the acquisition number provided, by 
#  first finding the closest scan number (scan index) and the finds the associated 
#  retention time.
########################################################################################
acqNumToRt <- function (xcmsRaw, aqNum) {
    xcmsRaw@scantime[as.integer(acqNumToScNum(xcmsRaw, aqNum))]
}

########################################################################################
# Input: xcmsRaw object, acquisition number
# Output: Scan number closest to acquisition number
# 
# This method finds the scan number (scan index) closest to the provided acquisition 
#  number by going through the list of acquisition numbers with associated scan numbers
#  and returning that scan number.
########################################################################################
acqNumToScNum <- function (xcmsRaw, aqNum) {
    which.min(abs(xcmsRaw@acquisitionNum - as.integer(aqNum)))
}

########################################################################################
# Input: mzXML path
# Output: xcmsRaw object
# 
# This method converts the mzXML at the path provided to a xcmsRaw object with 
#  "profstep=0".
########################################################################################
convert <- function (mzXMLPath) {
    return(xcmsRaw(mzXMLPath, profstep=0))
}

########################################################################################
# Input: xcmsRaw object, mass range (vector), RT range (vector)
# Output: EIC (raw)
# 
# This method geneerates an EIC with specified mass and retention time windows, and uses
#  the rawEIC method (not getEIC). Adjustments are made to the output from rawEIC to 
#  include retention times instead of scan numbers, by using the indices of scan numbers
#  as the indices of retention times (scantime).
# Final output that is returned has retention times and intensities associated with
#  them.
########################################################################################
getChromatogram <- function(xcmsRaw, mzrange, rtrange) {
    EIC <- rawEIC(xcmsRaw, mzrange, rtrange)
    EIC$rt <- unlist(lapply(EIC$scan, function(scan) xcmsRaw@scantime[scan]))
    EIC$scan <- NULL
    return(EIC)
}

########################################################################################
# Input: xcmsRaw object, mass list (vector), retention time, mass window (ppm), 
#  svm (for quality peak metrics), metric (boolean for calculating quality peak 
#  metrics), peakCalcLoader (boolean for only returning detected peaks, only used for
#  peakCalcLoader), rtw (default to constant at top of file, used for manual peak
#  detector)
# Output: List of peaks detected, standard deviation of peak apex retention times, and 
#  standard deviation of peak shapes
# 
# This follows what the "Run" method does in PEAKCALC.bas in visual basic.
# This method first runs selectMass on the masses provided, along with other parameters.
#  Then it goes through the peaks returned, either detected with centWave, or calculated
#  with getPeaks. Then it tries to use the retention time of a detected peak or the 
#  peak with the highest intensity, and redo the peak detection/calculation for the 
#  other peaks supplied. Once that is done, then it will calculate the quality peak 
#  metrics for those all peaks if required. Then the standard deviations are calcualted
#  with the retention times of the peaks (apex) and peak shape. This part is based on
#  the visual basic PEAKCALC.bas.
########################################################################################
peakCalc <- function(xcmsRaw, masslist, rt, massWindow, svmObject=NULL, metric=FALSE, peakCalcLoader=FALSE, rtw=RTWINDOW) {
    p <- lapply(masslist, function(mass) selectMass(xcmsRaw, mass, rt, massWindow, peakCalcLoader, rtw))
    #print(p)
	nPeaks <- length(p)
    peakAreas <- unlist(lapply(p, '[[', "peakArea"))

    if (!all(peakAreas > 0) || any(is.na(peakAreas))) {
        if (any(na.omit(peakAreas) > 0)) {
            maxPos <- which.max(peakAreas)
        } else {
            maxPos <- which.min(peakAreas)
            maxPos <- ifelse(length(maxPos) > 0, maxPos, -1)
        }

        if (maxPos > 0) {
            for (j in (1:nPeaks)[-maxPos]) {
                if (p[[j]][["peakArea"]] < 0 || is.na(p[[j]][["peakArea"]])) {
                    p[[j]] <- selectMass(xcmsRaw, masslist[j], rt, massWindow, peakCalcLoader, rtw)
                }
            }
        }
    }

    if (metric) {
        for (l in 1:nPeaks) {
            if (!is.na(p[[l]][["peakArea"]])) {
                #jmb added 150616 check if mass window was doubled to find peak, and adjust it if it was
			    if(p[[l]][["MZx2"]] == TRUE || p[[l]][["MZx2"]] == 1){
				    massWindow <- massWindow * 2
				}
				m <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], p[[l]][["peakLeftRT"]],p[[l]][["peakLeftsc"]], p[[l]][["peakRightRT"]], p[[l]][["peakRightsc"]], masslist[l], massWindow, p[[l]][["sigma"]])
#                if (length(m) < 1) {
                if (is.null(m[["metric"]])){ 
                    min <- head(xcmsRaw@scantime, n=1)
                    max <- tail(xcmsRaw@scantime, n=1)
                    peakLeft <- p[[l]][["peakApexRT"]] - 30
                    if (peakLeft < min) {
                        peakLeft <- min
                    }
					peakLeftsc = p[[l]][["peakLeftsc"]]
                    peakRight <- p[[l]][["peakApexRT"]] + 30
                    if (peakRight > max) {
                        peakRight <- max
                    }
					peakRightsc = p[[l]][["peakRightsc"]]
#                   tryCatch(p[[l]][["metric"]] <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], peakLeft, peakRight, masslist[l], massWindow), 
#                       error=function(ex) p[[l]][["metric"]] <- "no score")
                    m2 <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], peakLeft, peakLeftsc, peakRight, peakRightsc, masslist[l], massWindow, p[[l]][["sigma"]])
                    #use line below for calling metrics
					#m2 <- metrics(xcmsRaw, svmObject, rt, p[[l]][["peakApexRT"]], peakLeft,p[[l]][["peakLeftsc"]], peakRight, p[[l]][["peakRightsc"]], masslist[l], massWindow, p[[l]][["sigma"]])
                    if (!is.null(m2[["metric"]])) {
                        p[[l]][["metric"]] <- m2
                    } else {
                        p[[l]][["metric"]][["metric"]] <- "no score"
                    }
                } else {
                    p[[l]][["metric"]] <- m
                }
            } else {
                p[[l]][["metric"]][["metric"]] <- "no score"
            }
        }
    }
    stDevPeakApex <- standardDeviation(unlist(lapply(p, '[[', "peakApexRT")))

    lr <- rep(NA, nPeaks)
    for (k in 1:nPeaks) {
        if (is.na(p[[k]][["peakApexRT"]]) || p[[k]][["peakLeftRT"]] == p[[k]][["peakRightRT"]]) {
            lr[k] <- 0
        } else {
            lr[k] = (2 * p[[k]][["peakApexRT"]] - p[[k]][["peakLeftRT"]] - p[[k]][["peakRightRT"]]) / 
            (p[[k]][["peakRightRT"]] - p[[k]][["peakLeftRT"]])
        }
    }

    stDevPeakShape <- standardDeviation(lr)

    return(list(peaks=p, stDevPeakApex=stDevPeakApex, stDevPeakShape=stDevPeakShape))
}

########################################################################################
# Input: xcmsRaw object, mass, retention time, mass window (ppm), peakCalcLoader 
#  (boolean for only returning detected peaks, only used for peakCalcLoader), rtw 
#  (default to constant in top of file, used for manual peak detector)
# Output: Peak and associated information.
# 
# This follows what the "SelectMass" method does in PEAKCALC.bas in visual basic, but 
#  with modifications to run xcms and its peak detection, centWave.
# First this creates a mass range and a scan number range. These are used for the ROI
#  that will be used for centWave. Then centWave uses that roi, along with the 
#  parameters at the top of the script (and loops through multiple peak widths). The 
#  first run of centWave is with the whole scan range of the xcmsRaw object, while 
#  the second run of centWave is with the scan range calculated, and the third run is 
#  on double of that scan range. Then three sets of those can run, depending if the 
#  previous return a result or not. The first is the default centWave, the second is
#  with doubling the mass window, and checking if the peak detected fall within the
#  original mass window, and the third is with manually generating a ROI but with 
#  specific parameters (the main difference should be MINCENTROIDS, which in centWave
#  is defaulted to 4 as minimum). Then the outputs are merged, with the peaks detected 
#  with a larger peak width parameter at the top of the list. With that list, then it 
#  is filtered by retention times (peaks at the top of the list, or those that are 
#  detected with a larger peak width parameter), and NA's are removed. 
#  Now it goes through the case where centWave either detects something or not. If it 
#  does, then it will go through the list of peaks and find the one with the largest
#  intensity. If not, then it will run getPeaks with the center being the highest 
#  intensity point, and try to calculate a baseline adjusted intensity and signal to 
#  noise with a part of centWave method extracted from xcms source. 
# 
# http://metabolomics-forum.com/viewtopic.php?f=25&t=154&sid=1af1b9cff6303e147ab8e0c37fcd8626
########################################################################################
selectMass <- function(xcmsRaw, mass, rt, massWindow, peakCalcLoader=FALSE, rtw) {
    emptyMat <- matrix(c(rep(NA, 10)), nrow=1)
    noPeak <- list(peakArea=NA, peakApexRT=NA, peakLeftRT=NA, peakRightRT=NA, peakWidth=NA, sn=NA, metric=NA, MZx2 = FALSE)
    mzrange <- pmMz(mass, massWindow)
    rtrange <- pmRt(rt, rtw)
    scrange <- rtRangeToScanRange(xcmsRaw, rtrange)
	
    if (any(is.na(scrange))) return(noPeak)
    rtrange2 <- pmRt(rt, 2*rtw)
    scrange2 <- rtRangeToScanRange(xcmsRaw, rtrange2)
    roi <- list(list(mzmin=mzrange[1], mzmax=mzrange[2], scmin=scrange[1], scmax=scrange[2]))

    #####OLD############################################################################
    peaks <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roi, scrange, scrange2)
    filteredPeaks <- tryCatch(peaks[!duplicated(peaks[c("V4")]),], error = function(ex) NA)
	
    if (nrow(filteredPeaks) == 0 || is.na(filteredPeaks)) { # Second pass with double ppm and checking if any masses of peaks found fall within original mass window
		mzrange2 <- pmMz(mass, 2*massWindow)
        roi2 <- list(list(mzmin=mzrange2[1], mzmax=mzrange2[2], scmin=scrange[1], scmax=scrange[2]))

        peaks <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roi2, scrange, scrange2)
        if (nrow(peaks) > 0) {
            match1 <- intersect(which(peaks[,"V1"] >= mzrange[1]), which(peaks[,"V1"] <= mzrange[2])) # if mz of is within mzrange of original ppm
            match2 <- intersect(which(peaks[,"V2"] >= mzrange[1]), which(peaks[,"V2"] <= mzrange[2])) # if mzmin of is within mzrange of original ppm
            match3 <- intersect(which(peaks[,"V3"] >= mzrange[1]), which(peaks[,"V3"] <= mzrange[2])) # if mzmax of is within mzrange of original ppm
            indices <- unique(c(match1, match2, match3))
            # indices <- intersect(which(peaks[,"V2"] >= mzrange[1]), which(peaks[,"V3"] <= mzrange[2])) # if whole mz range of peak is within mzrange
            newPeaks <- peaks[indices,]
            filteredPeaks <- tryCatch(newPeaks[!duplicated(newPeaks[c("V4")]),], error = function(ex) NA)
            if(nrow(filteredPeaks) != 0 ){
             	filteredPeaks[["V23"]] <- TRUE #jmb added 150615--boolean signifies that mz was doubled
			}
		}
    } else {
	     filteredPeaks[["V23"]] <- FALSE #jmb added 150615
	}

	
    if (nrow(filteredPeaks) == 0 || is.na(filteredPeaks)) { # Third pass, manually detecting roi's and filtering out list
        sink(tempfile())
        # roisc <- tryCatch(xcms:::findmzROI(xcmsRaw, scanrange=scrange, dev=DEV, minCentroids=MINCENTROIDS, prefilter=PREFILTER), error = function(ex) matrix(nrow=0, ncol=0))
        roisc1 <- tryCatch(xcms:::findmzROI(xcmsRaw, scanrange=scrange, dev=DEV, minCentroids=MINCENTROIDS, prefilter=PREFILTER), error = function(ex) matrix(nrow=0, ncol=0))
        roisc2 <- tryCatch(xcms:::findmzROI(xcmsRaw, scanrange=scrange2, dev=DEV, minCentroids=MINCENTROIDS, prefilter=PREFILTER), error = function(ex) matrix(nrow=0, ncol=0))
        roisc <- unique(c(roisc1, roisc2))
        sink()
        if (length(roisc) > 0) {
            roidf <- data.frame(matrix(unlist(roisc), nrow=length(roisc), byrow=T))
            roiMatch1 <- intersect(which(roidf[,"X1"] >= mzrange[1]), which(roidf[,"X1"] <= mzrange[2]))
            roiMatch2 <- intersect(which(roidf[,"X2"] >= mzrange[1]), which(roidf[,"X2"] <= mzrange[2]))
            roiMatch3 <- intersect(which(roidf[,"X3"] >= mzrange[1]), which(roidf[,"X3"] <= mzrange[2]))
            roiIndices <- unique(c(roiMatch1, roiMatch2, roiMatch3))
            if (length(roiIndices) > 0) {
                newdf <- roidf[roiIndices,]
                roiL <- lapply(1:nrow(newdf), function(i) list(mzmin=newdf[i,"X2"], mzmax=newdf[i,"X3"], scmin=newdf[i,"X4"], scmax=newdf[i,"X5"]))
                peaks <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roiL, scrange, scrange2)  
			   if (nrow(peaks) > 0) {
                    match1 <- intersect(which(peaks[,"V1"] >= mzrange[1]), which(peaks[,"V1"] <= mzrange[2])) # if mz of is within mzrange of original ppm
                    match2 <- intersect(which(peaks[,"V2"] >= mzrange[1]), which(peaks[,"V2"] <= mzrange[2])) # if mzmin of is within mzrange of original ppm
                    match3 <- intersect(which(peaks[,"V3"] >= mzrange[1]), which(peaks[,"V3"] <= mzrange[2])) # if mzmax of is within mzrange of original ppm
                    # indices <- unique(c(match1, match2, match3))
                    indicesmz <- unique(c(match1, match2, match3))
                    match4 <- intersect(which(peaks[,"V4"] >= rtrange[1]), which(peaks[,"V4"] <= rtrange[2])) # if rt of is within rtrange of original rtw
                    match5 <- intersect(which(peaks[,"V5"] >= rtrange[1]), which(peaks[,"V5"] <= rtrange[2])) # if rtmin of is within rtrange of original rtw
                    match6 <- intersect(which(peaks[,"V6"] >= rtrange[1]), which(peaks[,"V6"] <= rtrange[2])) # if rtmax of is within rtrange of original rtw
                    indicesrt <- unique(c(match4, match5, match6))
                    indices <- intersect(indicesmz,indicesrt)
                    newPeaks <- peaks[indices,]
                    filteredPeaks <- tryCatch(newPeaks[!duplicated(newPeaks[c("V4")]),], error = function(ex) NA)
					if(nrow(filteredPeaks) != 0){
             	        filteredPeaks[["V23"]] <- FALSE #jmb added 150615--boolean signifies that mz was doubled	
			        }
                }
            }
        }
    }
    ####################################################################################

    #####ALL############################################################################
    # peaks1 <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roi, scrange, scrange2)
    # filteredPeaks1 <- tryCatch(peaks1[!duplicated(peaks1[c("V4")]),], error = function(ex) emptyMat)

    # mzrange2 <- pmMz(mass, 2*massWindow)
    # roi2 <- list(list(mzmin=mzrange2[1], mzmax=mzrange2[2], scmin=scrange[1], scmax=scrange[2]))

    # peaks2 <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roi2, scrange, scrange2)
    # filteredPeaks2 <- emptyMat
    # if (nrow(peaks2) > 0) {
    #     mzMatch11 <- intersect(which(peaks2[,"V1"] >= mzrange[1]), which(peaks2[,"V1"] <= mzrange[2])) # if mz of is within mzrange of original ppm
    #     mzMatch21 <- intersect(which(peaks2[,"V2"] >= mzrange[1]), which(peaks2[,"V2"] <= mzrange[2])) # if mzmin of is within mzrange of original ppm
    #     mzMatch31 <- intersect(which(peaks2[,"V3"] >= mzrange[1]), which(peaks2[,"V3"] <= mzrange[2])) # if mzmax of is within mzrange of original ppm
    #     indices1 <- unique(c(mzMatch11, mzMatch21, mzMatch31))
    #     # indices <- intersect(which(peaks[,"V2"] >= mzrange[1]), which(peaks[,"V3"] <= mzrange[2])) # if whole mz range of peak is within mzrange
    #     newPeaks1 <- peaks2[indices1,]
    #     filteredPeaks2 <- tryCatch(newPeaks1[!duplicated(newPeaks1[c("V4")]),], error = function(ex) emptyMat)
    # }

    # peaks3 <- NA
    # filteredPeaks3 <- emptyMat
    # sink(tempfile())
    # roisc <- tryCatch(xcms:::findmzROI(xcmsRaw, scanrange=scrange, dev=DEV, minCentroids=MINCENTROIDS, prefilter=PREFILTER), error = function(ex) matrix(nrow=0, ncol=0))
    # sink()
    # if (length(roisc) > 0) {
    #     roidf <- data.frame(matrix(unlist(roisc), nrow=length(roisc), byrow=T))
    #     roiMatch1 <- intersect(which(roidf[,"X1"] >= mzrange[1]), which(roidf[,"X1"] <= mzrange[2]))
    #     roiMatch2 <- intersect(which(roidf[,"X2"] >= mzrange[1]), which(roidf[,"X2"] <= mzrange[2]))
    #     roiMatch3 <- intersect(which(roidf[,"X3"] >= mzrange[1]), which(roidf[,"X3"] <= mzrange[2]))
    #     roiIndices <- unique(c(roiMatch1, roiMatch2, roiMatch3))
    #     if (length(roiIndices) > 0) {
    #         newdf <- roidf[roiIndices,]
    #         roiL <- lapply(1:nrow(newdf), function(i) list(mzmin=newdf[i,"X2"], mzmax=newdf[i,"X3"], scmin=newdf[i,"X4"], scmax=newdf[i,"X5"]))
    #         peaks3 <- cwMerge(xcmsRaw, PEAKWIDTHL, CWPPM, PREFILTER, SNTHRESH, INTEGRATE, MZDIFF, roiL, scrange, scrange2)
    #         if (nrow(peaks3) > 0) {
    #             mzMatch21 <- intersect(which(peaks3[,"V1"] >= mzrange[1]), which(peaks3[,"V1"] <= mzrange[2])) # if mz of is within mzrange of original ppm
    #             mzMatch22 <- intersect(which(peaks3[,"V2"] >= mzrange[1]), which(peaks3[,"V2"] <= mzrange[2])) # if mzmin of is within mzrange of original ppm
    #             mzMatch23 <- intersect(which(peaks3[,"V3"] >= mzrange[1]), which(peaks3[,"V3"] <= mzrange[2])) # if mzmax of is within mzrange of original ppm
    #             indices2 <- unique(c(mzMatch21, mzMatch22, mzMatch23))
    #             # indices <- intersect(which(peaks[,"V2"] >= mzrange[1]), which(peaks[,"V3"] <= mzrange[2])) # if whole mz range of peak is within mzrange
    #             newPeaks2 <- peaks3[indices2,]
    #             filteredPeaks3 <- tryCatch(newPeaks2[!duplicated(newPeaks2[c("V4")]),], error = function(ex) emptyMat)
    #         }
    #     }
    # }

    # filteredPeaksTemp <- rbind(filteredPeaks1, filteredPeaks2, filteredPeaks3)
    # filteredPeaks <- filteredPeaksTemp[rowSums(is.na(filteredPeaksTemp))!=10,]
    ####################################################################################

    if (nrow(filteredPeaks) == 0 || is.na(filteredPeaks)) { # if "No Peaks" or "omz == 0" (mass outside of raw file mass range) error
        if (peakCalcLoader) { # do not want pseudo peaks in detecting peaks with peakCalcLoader
            return(noPeak)
        }
        eic <- tryCatch(rawEIC(xcmsRaw, mzrange=mzrange, scanrange=scrange), error = function(ex) NA)
        if (all(eic$intensity == 0) || is.na(eic)) {
            return(noPeak)
        }
        maxIntIndex <- which.max(eic$intensity)
        maxInt <- eic$intensity[maxIntIndex]
        maxScan <- eic$scan[maxIntIndex]
        maxRt <- xcmsRaw@scantime[maxScan]
        maxRtrange <- pmRt(maxRt, rtw)
        scMaxrange <- rtRangeToScanRange(xcmsRaw, maxRtrange)
        roiMax <- list(list(mzmin=mzrange[1], mzmax=mzrange[2], scmin=scMaxrange[1], scmax=scMaxrange[2]))
        # baseline <- pseudoPeakCalcBaseline(xcmsRaw, ROI.list=roi, maxIntensity=maxInt)
        baselinesd <- pseudoPeakCalcBaseline2(xcmsRaw, ROI.list=roiMax, maxIntensity=maxInt)
        baseline <- baselinesd[["baseline"]]
        maxRtrange2 <- pmRt(maxRt, 2*rtw)
        eic2 <- rawEIC(xcmsRaw, mzrange=mzrange, rtrange=maxRtrange2)

        # pseudoPeak <- getPeaks(xcmsRaw, peakrange=matrix(c(mzmin=mzrange[1], mzmax=mzrange[2], rtmin=peakLeftRT, rtmax=peakRightRT), nrow=1), step=GSTEP)

            # pseudoPeak <- getPeaks(xcmsRaw, peakrange=matrix(c(mzmin=mzrange[1], mzmax=mzrange[2], rtmin=rtrange[1], rtmax=rtrange[2]), nrow=1))
        
        # if (pseudoPeak[1,][["into"]] == 0 || is.na(pseudoPeak)) {
        #     return(noPeak)
        # }

        # pseudoCalc <- pseudoPeakCalc(xcmsRaw, ROI.list=roiMax, maxIntensity=pseudoPeak[1,][["maxo"]], rawIntensity=pseudoPeak[1,][["into"]])

        # pseudoCalc <- pseudoPeakCalc2(xcmsRaw, ROI.list=roiMax, scanrange=scMaxrange, maxIntensity=maxInt, leftsc=eic2$scan[peakLeft], rightsc=eic2$scan[peakRight])
        # pseudoCalc <- pseudoPeakCalc2(xcmsRaw, ROI.list=roiMax, maxIntensity=maxInt, leftsc=eic2$scan[peakLeft], rightsc=eic2$scan[peakRight])
        # print(c(peakLeftsc,peakRightsc))

        tryCatch({

        maxIndex <- which(eic2$intensity == maxInt)
        peakLeft <- maxIndex - 1
        underBaseCount <- 0

        for (i in (maxIndex-1):1) {
            if (is.na(eic2$intensity[i])) {
                break
            }
            peakLeft <- i
            if (eic2$intensity[i] < baseline) {
                # break
                underBaseCount <- underBaseCount + 1
            } else {
                underBaseCount <- 0
            }
            if (underBaseCount >= BASELINECOUNT) {
                break
            }
        }
        peakRight <- maxIndex + 1
        underBaseCount <- 0

        for (i in (maxIndex+1):length(eic2$intensity)) {
            if (is.na(eic2$intensity[i])) {
                break
            }
            peakRight <- i
            if (eic2$intensity[i] < baseline) {
                # break
                underBaseCount <- underBaseCount + 1
            } else {
                underBaseCount <- 0
            }
            if (underBaseCount >= BASELINECOUNT) {
                break
            }
        }
        peakLeftsc <- eic2$scan[peakLeft]
        peakLeftRT <- xcmsRaw@scantime[peakLeftsc]
        peakRightsc <- eic2$scan[peakRight]
        peakRightRT <- xcmsRaw@scantime[peakRightsc]

            pseudoCalc <<- pseudoPeakCalc2(xcmsRaw, mzrange=mzrange, scanrange=c(peakLeftsc,peakRightsc), maxIntensity=maxInt, baseline=baseline, sdnoise=baselinesd[["sdnoise"]])
        }, error = function(e) {
            pseudoCalc <<- NA
        })

        suppressWarnings(if (is.na(pseudoCalc)) {
            return(noPeak)
        })

        #pseudoCalc <- pseudoPeakCalc2(xcmsRaw, mzrange=mzrange, scanrange=c(peakLeftsc,peakRightsc), maxIntensity=maxInt, baseline=baseline, sdnoise=baselinesd[["sdnoise"]])
        
        return(list(peakArea=-pseudoCalc[["intb"]], peakApexRT=maxRt, peakLeftRT=peakLeftRT, peakLeftsc = peakLeftsc,
                    peakRightRT=peakRightRT, peakRightsc = peakRightsc, peakWidth=peakRightRT-peakLeftRT, 
                    sn=pseudoCalc[["sn"]], metric=NA, sigma=NA, MZx2=FALSE)) #jmb this needs to be fixed--sigma shouldn't be NA
    } else {
        # pl <- unlist(filteredPeaks[which.min(abs(filteredPeaks[,"V4"] - as.numeric(rt))),]) # Pick peak with rt closest to input rt
		pl <- unlist(filteredPeaks[which.max(filteredPeaks[,"V8"]),]) # Pick peak with highest intensity
		return(list(peakArea=pl[["V8"]], peakApexRT=pl[["V4"]], peakLeftRT=pl[["V5"]], peakLeftsc=pl[["V19"]],
		  peakRightRT=pl[["V6"]], peakRightsc=pl[["V20"]], peakWidth=pl[["V6"]]-pl[["V5"]], sn=pl[["V10"]],
		  metric=NA, sigma=pl[["V13"]], MZx2=pl[["V23"]]))
    }
}

########################################################################################
# Input: xcmsRaw object, peak width , mass window, prefilter, snthresh, integrate, 
#  mzdiff, roi, scan range, scan range 2 (doubled of previous)
# Output: peaks as matrix (centwave output)
# 
# This method will run three sets of centwave, as stated in selectMass. The first one 
#  is with the normal ROI and global variables, the second is with a scan range, and
#  the third one is with doubled the original scanrange. A matrix is then returned,
#  of detected peaks.
########################################################################################
cwMerge <- function(xcmsRaw, peakwidths, ppm, prefilter, snthresh, integrate, mzdiff, roi, scrange, scrange2) {
    emptyMat <- matrix(c(rep(NA, 22)), nrow=1)
	sink(tempfile())
    cw <- lapply(peakwidths, function(w) {
        mat <- rbind(tryCatch(suppressWarnings(findPeaks.centWave(xcmsRaw, ppm=ppm, peakwidth=w, prefilter=prefilter, fitgauss=TRUE,
                                                        snthresh=snthresh, integrate=integrate, mzdiff=mzdiff, 
                                                        ROI.list=roi, verbose.columns=TRUE)), error = function(ex) {emptyMat}), 
        tryCatch(suppressWarnings(findPeaks.centWave(xcmsRaw, ppm=CWPPM, peakwidth=w, prefilter=prefilter, fitgauss=TRUE,
                                                        snthresh=snthresh, integrate=integrate, mzdiff=mzdiff, 
                                                        ROI.list=roi, scanrange=scrange, verbose.columns=TRUE)), error = function(ex) {emptyMat}),
        tryCatch(suppressWarnings(findPeaks.centWave(xcmsRaw, ppm=CWPPM, peakwidth=w, prefilter=prefilter, fitgauss=TRUE,
                                                        snthresh=snthresh, integrate=integrate, mzdiff=mzdiff, 
                                                        ROI.list=roi, scanrange=scrange2, verbose.columns=TRUE)), error = function(ex) {emptyMat}))
        matc <- matrix(mat[rowSums(is.na(mat))!=22,], ncol=22)
        if (nrow(matc) == 0) {
            return(NA)
        } else {
            return(matc)
        }
        })
    sink()
    peaks <- as.data.frame(do.call(rbind, rev(cw[!is.na(cw)])), ncol=22)
    return(peaks)
}

########################################################################################
# Input: xcmsRaw object, peak width (global variable), signal to noise threshold 
#  (global), prefilter (global), scanrange (NULL), ROI.list (NULL), maxIntensity (from 
#  getPeaks, maximum intensity of given range (EIC)), rawIntensity (into from getPeaks)
# Output: Signal to noise, intb (baseline adjusted intensity)
# 
# This method is based on the source for xcms's centWave peak detection method. Signal
#  to noise is calculated in the same manner as centWave would. Peak area (intensity) 
#  had to be adjusted as a peak was not detected. The major adjustment made was there
#  is no narrowing down of the peak width, but instead uses the whole scan number range
#  provided in the ROI list. 
# Based on code from XCMS 1.40.0
########################################################################################
pseudoPeakCalc <- function(object, peakwidth=PEAKWIDTHG, snthresh=SNTHRESH, prefilter=PREFILTER, scanrange=numeric(), ROI.list=list(), maxIntensity, rawIntensity) {
    maxint <- unname(maxIntensity)

    # Code from centwave, modified to calculate "sn" and "intb" from getPeaks() output
    scanrange.old <- scanrange
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    # if (!(identical(scanrange.old,scanrange)) && (length(scanrange.old) >0))
    #     cat("Warning: scanrange was adjusted to ",scanrange,"\n")

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(object@scantime))) / 2)

    if (length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]

    # if (length(scalerange) < 1)
    #     stop("No scales ? Please check peak width!\n")

    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
    scales <- scalerange;

    minPeakWidth <-  scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2);

    peaklist <- list()
    scantime <- object@scantime
    Nscantime <- length(scantime)

    feat <- ROI.list[[1]]
    N <- feat$scmax - feat$scmin + 1

    peaks <- peakinfo <- NULL

    mzrange <- c(feat$mzmin,feat$mzmax)
    sccenter <- feat$scmin[1] + floor(N/2) - 1
    scrange <- c(feat$scmin,feat$scmax)
    ## scrange + noiserange, used for baseline detection and wavelet analysis
    sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
    eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
    d <- eic$intensity
    td <- sr[1]:sr[2]
    scan.range <- c(sr[1],sr[2])
    ## original mzROI range
    mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
    omz <- xcms:::rawMZ(object,mzrange=mzrange,scanrange=scrange)
    # if (all(omz == 0))
    #     stop("centWave: debug me: (omz == 0)?\n")
    od  <- mzROI.EIC$intensity
    otd <- mzROI.EIC$scan
    # if (all(od == 0))
    #     stop("centWave: debug me: (all(od == 0))?\n")

    ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
    ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
    fd <- d[match(ftd,td)]

    ## 1st type of baseline: statistic approach
    if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
        noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity else
    noised <- d;
    ## 90% trimmed mean as first baseline guess
    noise <- xcms:::estimateChromNoise(noised, trim=0.05, minPts=3*minPeakWidth)

    # ## any continuous data above 1st baseline ?
    # if (!continuousPtsAboveThreshold(fd,threshold=noise,num=minPtsAboveBaseLine))
    #     next;

    ## 2nd baseline estimate using not-peak-range
    lnoise <- xcms:::getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime, threshold=noise,num=minPtsAboveBaseLine)

    ## Final baseline & Noise estimate
    baseline <- max(1,min(lnoise[1],noise))
    sdnoise <- max(1,lnoise[2])
    sdthr <-  sdnoise * snthresh

    sn <- round((maxint - baseline) / sdnoise)

    lm <- c(1, length(d))

    into <- sum(d[lm[1]:lm[2]])
    ratio <- unname(rawIntensity / into)
    newpwid <- ifelse(ratio > 1, ratio, 1)

    db <-  d[lm[1]:lm[2]] - baseline
    intb <- newpwid*sum(db[db>0])

    return(list(intb=intb, sn=sn))
}

########################################################################################
# Input: xcmsRaw object, peak width (global variable), signal to noise threshold 
#  (global), prefilter (global), scanrange (NULL), ROI.list (NULL), maxIntensity (from 
#  getPeaks, maximum intensity of given range (EIC))
# Output: baseline
# 
# This method is based on the source for xcms's centWave peak detection method, and is 
#  the same as the above method but does not continue on after a base line is calcualted
#  and returned.
# Based on code from XCMS 1.40.0
########################################################################################
pseudoPeakCalcBaseline <- function(object, peakwidth=PEAKWIDTHG, snthresh=SNTHRESH, prefilter=PREFILTER, scanrange=numeric(), ROI.list=list(), maxIntensity) {
    maxint <- unname(maxIntensity)

    # Code from centwave, modified to calculate baseline
    scanrange.old <- scanrange
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(object@scantime))) / 2)

    if (length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]

    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
    scales <- scalerange;

    minPeakWidth <-  scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2);

    peaklist <- list()
    scantime <- object@scantime
    Nscantime <- length(scantime)

    feat <- ROI.list[[1]]
    N <- feat$scmax - feat$scmin + 1

    peaks <- peakinfo <- NULL

    mzrange <- c(feat$mzmin,feat$mzmax)
    sccenter <- feat$scmin[1] + floor(N/2) - 1
    scrange <- c(feat$scmin,feat$scmax)
    ## scrange + noiserange, used for baseline detection and wavelet analysis
    sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
    eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
    d <- eic$intensity
    td <- sr[1]:sr[2]
    scan.range <- c(sr[1],sr[2])
    ## original mzROI range
    mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
    omz <- xcms:::rawMZ(object,mzrange=mzrange,scanrange=scrange)
    od  <- mzROI.EIC$intensity
    otd <- mzROI.EIC$scan

    ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
    ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
    fd <- d[match(ftd,td)]

    ## 1st type of baseline: statistic approach
    if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
        noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity else
    noised <- d;
    ## 90% trimmed mean as first baseline guess
    noise <- xcms:::estimateChromNoise(noised, trim=0.05, minPts=3*minPeakWidth)

    ## 2nd baseline estimate using not-peak-range
    lnoise <- xcms:::getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime, threshold=noise,num=minPtsAboveBaseLine)

    ## Final baseline & Noise estimate
    baseline <- max(1,min(lnoise[1],noise))
    sdnoise <- max(1,lnoise[2])
    sdthr <-  sdnoise * snthresh

    return(baseline)
}

pseudoPeakCalcBaseline2 <- function(object, peakwidth=PEAKWIDTHG, snthresh=SNTHRESH, prefilter=PREFILTER, scanrange=numeric(), ROI.list=list(), maxIntensity) {
    maxint <- unname(maxIntensity)

    # Code from centwave, modified to calculate baseline
    scanrange.old <- scanrange
    if (length(scanrange) < 2)
        scanrange <- c(1, length(object@scantime))
    else
        scanrange <- range(scanrange)

    scanrange[1] <- max(1,scanrange[1])
    scanrange[2] <- min(length(object@scantime),scanrange[2])

    ## Peak width: seconds to scales
    scalerange <- round((peakwidth / mean(diff(object@scantime))) / 2)

    if (length(z <- which(scalerange==0)))
        scalerange <- scalerange[-z]

    if (length(scalerange) > 1)
        scales <- seq(from=scalerange[1], to=scalerange[2], by=2)  else
    scales <- scalerange;

    minPeakWidth <-  scales[1];
    noiserange <- c(minPeakWidth*3, max(scales)*3);
    maxGaussOverlap <- 0.5;
    minPtsAboveBaseLine <- max(4,minPeakWidth-2);
    minCentroids <- minPtsAboveBaseLine ;
    scRangeTol <-  maxDescOutlier <- floor(minPeakWidth/2);

    peaklist <- list()
    scantime <- object@scantime
    Nscantime <- length(scantime)

    feat <- ROI.list[[1]]
    N <- feat$scmax - feat$scmin + 1

    peaks <- peakinfo <- NULL

    mzrange <- c(feat$mzmin,feat$mzmax)
    sccenter <- feat$scmin[1] + floor(N/2) - 1
    scrange <- c(feat$scmin,feat$scmax)
    ## scrange + noiserange, used for baseline detection and wavelet analysis
    sr <- c(max(scanrange[1],scrange[1] - max(noiserange)),min(scanrange[2],scrange[2] + max(noiserange)))
    eic <- rawEIC(object,mzrange=mzrange,scanrange=sr)
    d <- eic$intensity
    td <- sr[1]:sr[2]
    scan.range <- c(sr[1],sr[2])
    ## original mzROI range
    mzROI.EIC <- rawEIC(object,mzrange=mzrange,scanrange=scrange)
    omz <- xcms:::rawMZ(object,mzrange=mzrange,scanrange=scrange)
    od  <- mzROI.EIC$intensity
    otd <- mzROI.EIC$scan

    ##  scrange + scRangeTol, used for gauss fitting and continuous data above 1st baseline detection
    ftd <- max(td[1], scrange[1] - scRangeTol) : min(td[length(td)], scrange[2] + scRangeTol)
    fd <- d[match(ftd,td)]

    ## 1st type of baseline: statistic approach
    if (N >= 10*minPeakWidth)  ## in case of very long mass trace use full scan range for baseline detection
        noised <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity else
    noised <- d;
    ## 90% trimmed mean as first baseline guess
    noise <- xcms:::estimateChromNoise(noised, trim=0.05, minPts=3*minPeakWidth)

    ## 2nd baseline estimate using not-peak-range
    lnoise <- xcms:::getLocalNoiseEstimate(d,td,ftd,noiserange,Nscantime, threshold=noise,num=minPtsAboveBaseLine)

    ## Final baseline & Noise estimate
    baseline <- max(1,min(lnoise[1],noise))
    sdnoise <- max(1,lnoise[2])
    sdthr <-  sdnoise * snthresh

    return(list(baseline=baseline, sdnoise=sdnoise))
}

pseudoPeakCalc2 <- function(object, mzrange, scanrange, maxIntensity, baseline, sdnoise) {
    maxint <- unname(maxIntensity)

    sn <- round((maxint - baseline) / sdnoise)

    pwid <- (object@scantime[scanrange[2]] - object@scantime[scanrange[1]])/(scanrange[2] - scanrange[1])
    if (is.na(pwid))
        pwid <- 1

    tmpeic <- rawEIC(object,mzrange=mzrange,scanrange=scanrange)$intensity

    # into <- sum(d[lm[1]:lm[2]])
    # ratio <- unname(rawIntensity / into)
    # newpwid <- ifelse(ratio > 1, ratio, 1)

    # db <-  d[lm[1]:lm[2]] - baseline
    db <- tmpeic - baseline
    intb <- pwid*sum(db[db>0])

    return(list(intb=intb, sn=sn))
}

# findPeaks <- function(x, thresh=0) {
#     pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) < 0) + 2
#     if(!missing(thresh)) {
#      pks[x[pks-1]-x[pks] > thresh]
#     } else pks
# }

# findValleys <- function(x, thresh=0) {
#     pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) > 0) + 2
#     if(!missing(thresh)) {
#      pks[x[pks-1]-x[pks] > thresh]
#     } else pks
# }

########################################################################################
# Input: List
# Output: Standard deviation (number)
# 
# This method will calculate the standard deviation of a list (of size greater than 2)
#  in the same manner as the visual basic PEAKCALC.bas did. 
########################################################################################
standardDeviation <- function(list) {
    dblAnswer <- 0
    lengthList <- length(list[na.omit(list) > 0])

    if (lengthList > 1) {
        dblSum <- sum(list, na.rm=TRUE)
        dblMean <- dblSum / lengthList
        dblSumSqdDevs <- sum((list[list > 0] - dblMean) ^ 2, na.rm=TRUE)
        dblAnswer <- sqrt(dblSumSqdDevs / (lengthList - 1))
    }
    return(dblAnswer)
}


########################################################################################
# Metrics ##############################################################################
# ms2Call <- function (rTime, rt) {
#     return(abs(rt-mean(rTime)) / 60)
# }

ms2Call <- function (cwrt, rt) {
	#changed 11032016 to return value in seconds
    return(abs(rt-cwrt))
}

noiseDetector <- function (rTime, iList, rEdgev, lEdgev) {
    rightDetect <- tail(rTime, n=1)
    leftDetect <- head(rTime, n=1)
    maxiList <- max(iList)*0.2
    noiseScore <- 0
    index <- 1

    while (leftDetect < lEdgev) {
        if (iList[match(leftDetect, rTime)] > maxiList) {
            noiseScore = noiseScore + 1
        }
        index = index + 1
        leftDetect = rTime[index]
    }

    index <- length(rTime)

    while (rightDetect > rEdgev) {
        if (iList[match(rightDetect, rTime)] > maxiList) {
            noiseScore = noiseScore + 1
        }
        index = index - 1
        rightDetect = rTime[index]
    }
    return(noiseScore)
}

normalityTests <- function (iList, rEdge, lEdge) {
    tryCatch({
        return(shapiro.test(iList[lEdge:(rEdge+1)])$p.value)
    }, error = function(ex) {
        return(1)
    })
}

derivativeChanges <- function (rTime, iList, rEdge, lEdge) {
    curr.slope = ((iList[2])-(iList[1]))/((rTime[2])-(rTime[1]))
    num.changes = 0
    curr.intensity = iList[2]
    negthreshold = 0.5
    posthreshold=1.5

    tryCatch({
        for (i in (lEdge+2):rEdge-1) {
            new.slope = ((iList[i+1])-(iList[i]))/((rTime[i+1])-(rTime[i]))
            new.intensity = iList[i+1]
            if (new.slope * curr.slope < 0) {
                if (((new.slope < 0) & (new.intensity <= curr.intensity/negthreshold)) |
                            ((new.slope > 0) & (new.intensity >= curr.intensity*posthreshold))) {
                    num.changes = num.changes + 1
                    curr.intensity = new.intensity
                    curr.slope = new.slope
                }
            } else if (new.slope * curr.slope > 0) {
                if (((new.slope > 0) & (new.intensity > curr.intensity)) |
                            (( new.slope < 0) & (new.intensity < curr.intensity))) {
                    curr.intensity = new.intensity
                    curr.slope = new.slope
                }
            } else if (new.slope * curr.slope == 0) {
                if (new.intensity != curr.intensity) {
                    curr.intensity = new.intensity
                    num.changes = num.changes+1
                    curr.slope = new.slope
                }
            }
        }
        return(num.changes)
    }, error = function(ex) {
        return(0)
    })
}

zeroCounts <- function (iList, rEdge, lEdge) {
    num.changes = 0
    initial.intensity = iList[lEdge]

    tryCatch({
        for (i in lEdge:rEdge) {
            if ((initial.intensity == 0) & (num.changes == 0)) {
                if (iList[i+1]!= 0) {
                    num.changes = 1
                }
            } else {
                if ((iList[i] != 0) & (iList[i+1] == 0)) {
                    num.changes = num.changes +1
                }
            }
        }
        return(num.changes)
    }, error = function(ex) {
        return(0)
    })
}

maxIntensity <- function (iList, rEdge, lEdge) {
    return(max(iList[lEdge:(rEdge+1)], na.rm=TRUE))	#jmb 150507: added na.rm=TRUE
}

outsideIntensitySum <- function(iList, rEdge, lEdge, rEdgev, lEdgev){
    return(sum(append(iList[1:lEdge],iList[rEdge:length(iList)]))/rtWindow(rEdgev,lEdgev)) #2 * rtWindow size?
}

# insideIntensitySum <- function(iList, rEdge, lEdge, rEdgev, lEdgev){
    # return(sum(iList[lEdge,rEdge])/rtWindow(rEdgev,lEdgev))
# }

rtWindow <- function (rEdgev, lEdgev) {
	#changed 11032016 to return rtWindow in seconds
    return(rEdgev - lEdgev)
}

kurtosis <- function (iList){
    iUnList <- unlist(iList)
    m <- mean(iUnList, na.rm=TRUE)
    N <- length(iUnList)
    k <- N * sum((iUnList - m)^4)/(sum((iUnList - m)^2)^2)
    return(k)
}

#####Can deprecate this! skewness() is already implemented in e1071
skew <- function(iList){
    iUnList <- unlist(iList)
    m <- mean(iUnList, na.rm=TRUE)
    N <- length(iUnList)
    skewness <- (sum((iUnList - m)^3)/N)/(sum((iUnList - m)^2)/N)^(3/2)
    return(skewness)
}

DistToInflectionPoint <- function(iList, threshold, rEdge, lEdge, from = "L", to = "R"){
    intensity <- iList
    intensity <- intensity - threshold
    distance <- 0
    if(from == "L" & to == "R"){
        #Go through table and look for the first inflection point from L->R and R->L
        for (i in (lEdge):rEdge) {
          #If intensity is above the defined baseline and upcoming scans are also above it, then break
            if (intensity[i] > 0) {
                if(i+2 >= rEdge || intensity[i+2] >= 0){break}
        
            } else { distance <- distance + 1}
        }
    } else if (from == "R" & to == "L"){
		for (i in (rEdge):lEdge) {
		  #If intensity is above the defined baseline and upcoming scans are also above it, then break
			if (intensity[i] > 0) {
				if(i-2 <= lEdge || intensity[i-2] >= 0){break}
		
			} else { distance <- distance + 1}
		}
    } else {print("ERROR: invalid arguments provided")}
    #Convert number of zero/baseline scans to a percentage of the total rtwindow size before returning result
    return(distance/(rEdge-lEdge))
}

#Zhang and Zao, BMC Bioinformatics (2014); Choi, Biotechnology and Bioprocess Engineering (2004)
peakSharpness <- function(iList){
	intensity <- iList
	maxi <- which.max(intensity)
	lEdge <- 1
	rEdge <- length(intensity)
	rise <- 0
	if(maxi != 1){
		for(i in (lEdge+1):maxi){
			if(intensity[i-1] > 0){
				rise <- rise + ((intensity[i]-intensity[i-1])/intensity[i-1])
			} else {
				rise <- rise + intensity[i]
			}
		}
	}
	fall <- 0
	if(maxi != rEdge){
		for(i in maxi:(rEdge-1)){
			if(intensity[i+1] > 0){
				fall <- fall + ((intensity[i]-intensity[i+1])/intensity[i+1])
			} else {
				fall <- fall + intensity[i]
			}
		}
	}
	return(rise+fall)
}

#Zhang and Zao, BMC Bioinformatics (2014)
TPASR <- function(iList, rtL, rtR){
	intensity <- iList
	lEdge <- 1
	rEdge <- length(intensity)
	tpa <- 0.5 * (rtR-rtL) * max(intensity)
	rpa <- sum(intensity)
	tpasr <- abs(tpa-rpa)/tpa
	return(tpasr)
}

#Zhang and Zao, BMC Bioinformatics (2014)
peakSignificance <- function(iList,split){
	intensity <- iList
	if(length(intensity) < split){
		return(max(intensity)/min(intensity))
	}
	regionSize <- length(intensity)/split
	lWindow <- intensity[1:regionSize]
	midWindow <- intensity[(2*regionSize):(3*regionSize)]
	rWindow <- intensity[(4*regionSize):length(intensity)]

	midMean <- mean(midWindow)
	lrMean <- mean(c(lWindow, rWindow))
	if(lrMean < 1){
		lrMean <- 1
	}
	return(midMean/lrMean)
}

# FWHM <- function (rTime, iList, lEdge, rEdge) {
	# iList = iList[lEdge:rEdge]
	# rTime = rTime[lEdge:rEdge]
	# halfmax <- max(iList)/2
    # max_index <- which.max(iList)
    
	# left = max_index
	# lastdiff = Inf
	# i = max_index-1
	# while(i >= 1){
		# currdiff = abs(iList[i]-halfmax)
		# if(iList[i] > halfmax | currdiff < lastdiff){
			# left = i
			# lastdiff = currdiff
		# } else{
			# break
		# }
		# i = i-1
	# }

	# right = max_index
	# lastdiff = Inf
	# i = max_index+1
	# while(i <= length(iList)){
		# currdiff = abs(iList[i]-halfmax)
		# if(iList[i] > halfmax | currdiff < lastdiff){
			# right = i
			# lastdiff = currdiff
		# } else{
			# break
		# }
		# i = i+1
	# }

    # return(rTime[right]-rTime[left])
# }

FWHM = function(sigma, lrt, lsc, rrt, rsc){
	perc = sigma/(rsc-lsc)
	sigma = (rrt-lrt) * perc
	FWHM = 2 * sqrt(2 * log(2)) * sigma
	return(FWHM)
}

svmObjectm <- function (path) {
    # #dat <- read.csv(path, colClasses=c('NULL',NA,NA,NA,NA,'NULL',NA,NA,NA,NA))
    # dat <- read.table(path, colClasses=c('NULL',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),header=TRUE) #14 col

    # BScut <- function (x, cut) {
        # x <- ifelse(x <= cut, 1, 0)
    # }

    # cutoff <- 3
    # labels <- unlist(lapply(dat[,2], BScut, cutoff))

    # gamma <- 0.15
    # cost <- 0.55
	
    # return(svm(Rank ~ ., data = dat, cost = cost, gamma = gamma))
	
	dat <- read.table(path, colClasses=c('NULL','NULL',NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),header=TRUE, sep="\t") #15 col
	ranks <- as.factor(dat[,1])
	attribs <- dat[,2:ncol(dat)]
	
	gamma <- 0.15
	cost <- 0.55
	
	#weights are set to the inverse of that label's frequency in the data 
	nWeight <- sum(dat[,1]==1)/nrow(dat)
	pWeight <- sum(dat[,1]==0)/nrow(dat)
	
	return(svm(attribs, ranks, type="C-classification", cost=cost, gamma=gamma, class.weights=c("0"=nWeight,"1"=pWeight)))
}

metrics <- function (xcmsRaw, svmObject, rt, cwrt, rtmin, scmin, rtmax, scmax, mz, ppm, sigma) {
	EIC <- getChromatogram(xcmsRaw, mzrange=pmMz(mz, ppm), rtrange=c(rtmin-(rtmax-rtmin), rtmax+(rtmax-rtmin)))
    rTime <- EIC$rt
    iList <- EIC$intensity
    rEdge <- which.min(abs(rTime-rtmax))
    lEdge <- which.min(abs(rTime-rtmin))
	#If all intensities in SIC are 0, this likely means that mass window was doubled to find peak.
	#For now, allow these peaks, but calculate the quality metrics on the expanded SIC with wider MZ
	if(all(unlist(iList)==0)){
	    return(list(metric=NULL,ms2Call=NA,noise=NA,normality=NA,
        derivativeCount=NA,zeroCount=NA,maxIntensity=NA,
        rtWindow=NA,kurtosis=NA,skew=NA,peakSig=NA,
		peakSharp=NA,TPASR=NA,FWHM=NA
		))
	}

    if (rEdge == lEdge) {
        if (lEdge > 1) {
            lEdge <- lEdge - 1
        }
        if (rEdge < length(rTime)) {
            rEdge <- rEdge + 1
        }
    }
    rEdgev <- rTime[rEdge]
    lEdgev <- rTime[lEdge]

    ms2CallV <- ms2Call(cwrt, rt)
    noiseDetectorV <- noiseDetector(rTime, iList, rEdgev, lEdgev)
    normalityTestsV <- normalityTests(iList, rEdge, lEdge)
    derivativeChangesV <- derivativeChanges(rTime, iList, rEdge, lEdge)
    zeroCountsV <- zeroCounts(iList, rEdge, lEdge)
    maxIntensityV <- maxIntensity(iList, rEdge, lEdge)
    rtWindowV <- rtWindow(rtmax, rtmin)
    kurtosisV <- kurtosis(iList)
    skewV <- skew(iList)
    outsideIntensitySumV <- outsideIntensitySum(iList, rEdge, lEdge, rEdgev, lEdgev)
    distToInflL2RV <- DistToInflectionPoint(iList, 1000, rEdge, lEdge, from = "L", to = "R")
    distToInflR2LV <- DistToInflectionPoint(iList, 1000, rEdge, lEdge, from = "R", to = "L")
    peakSignificanceV <- peakSignificance(iList, 5)
    peakSharpnessV <- peakSharpness(iList)
    TPASRV <- TPASR(iList, lEdge, rEdge)
    tryCatch(FWHMV <- FWHM(sigma, lEdge, scmin, rEdge, scmax), 
      error=function(ex) FWHMV <- NA)
#    insideIntensitySumV <- insideIntensitySum(iList, rEdge, lEdge, rEdgev, lEdgev)
    # return(predict(svmObject, matrix(c(ms2Call(cwrt, rt), 
    #     noiseDetector(rTime, iList, rEdgev, lEdgev), 
    #     normalityTests(iList, rEdge, lEdge), 
    #     derivativeChanges(rTime, iList, rEdge, lEdge), 
    #     zeroCounts(iList, rEdge, lEdge), 
    #     maxIntensity(iList, rEdge, lEdge), 
    #     rtWindow(rEdgev, lEfdgev)), nrow=1, ncol=7), decision.values=FALSE))    

	metricV <- predict(svmObject, matrix(c(ms2CallV, 
        noiseDetectorV, 
        normalityTestsV, 
        derivativeChangesV, 
        zeroCountsV, 
        maxIntensityV, 
        rtWindowV,
		kurtosisV,
		skewV,
		outsideIntensitySumV,
		distToInflL2RV,
		distToInflR2LV
		), nrow=1, ncol=12), decision.values=FALSE)

    if (length(metricV) < 1) {
        metricV <- NULL
    }
    return(list(metric=as.character(metricV),ms2Call=ms2CallV,noise=noiseDetectorV,normality=normalityTestsV,
        derivativeCount=derivativeChangesV,zeroCount=zeroCountsV,maxIntensity=maxIntensityV,
        rtWindow=rtWindowV,kurtosis=kurtosisV,skew=skewV
		#,outsideIntensitySum=outsideIntensitySumV,distToInflL2RV,distToInflR2LV
		,peakSig=peakSignificanceV,peakSharp=peakSharpnessV,TPASR=TPASRV,FWHM=FWHMV
		))
}