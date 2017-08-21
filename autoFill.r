#!/usr/bin/Rscript
# @aly

########################################################################################
# To make this script run, one must provide three arguments in the following order:
#  1. Config file path (config file holds constants, like paths and file names)
#  2. Replicate number
#  3. SILAC or LABELFREE string, to tell this script to run what type of quantitation
#  4. PARALLEL or NOPARALLEL string, indicator for multiprocessing
# Also the following packages must be installed:
#  parallel (current versions of R should have it bundled)
#  XML (from CRAN)
#  xcms (from bioconductor)
#  e1071 (from CRAN)
#  nortest (from CRAN)
########################################################################################

suppressWarnings(suppressPackageStartupMessages({
    library(tcltk)
    library(parallel)
    library(XML)
    library(xcms)
    library(e1071)
    library(nortest)
    }))

#Suppress scientific notation
options(scipen=999)

########################################################################################
# Names of xcms methods script and training set csv for svm. Location of them should be in 
#  same folder as script.
########################################################################################
xcmsScript <<- "xcms2.r"
trainingSet <<- "SVMTrainingSet.txt"
progressFolder <<- "progress/"

########################################################################################
# Input: xml path
# Output: parsed xml as matrix
# 
# This method will preallocate and parse filemaker style xml's into a matrix. It first
#  looks for the "FIELD" tags and stores the lengths of the repeating fields. When it
#  encounters the "RESULTSET" tag, it will generate a list of lists and preallocate them
#  for filling in later on. In the attributes of "RESULTSET", it will use the number
#  (number of records) as the size of the outer list. Only entries in the "DATA" tags
#  are stored. At the end, the list of lists is flattened into a matrix.
# This uses SAX style xml parsing. 
########################################################################################
parseFMXML <- function(path) {
    fieldCounter <- 0
    rowCounter <- 0
    colCounter <- 0
    datCounter <- 0
    colBool <- FALSE
    datBool <- FALSE
    fieldSizes <- list()
    parsedXML <- list()
    xmlEventParse(
        path, 
        handlers = list(
            startElement = function(name, attr) {
                if (name == "FIELD") {
                    fieldCounter <<- fieldCounter + 1
                    fieldSizes[[fieldCounter]] <<- as.numeric(attr[[2]])
                }
                if (name == "RESULTSET") {
                    rowCounter <<- 0
                    parsedXML <<- lapply(1:attr, function(x) lapply(fieldSizes, function(y) rep(NA, y)))
                }
                if (name == "ROW") {
                    rowCounter <<- rowCounter + 1
                    colCounter <<- 0
                    datCounter <<- 0
                }
                if (name == "COL") {
                    colBool <<- TRUE
                    colCounter <<- colCounter + 1
                }
                if (name == "DATA") {
                    datBool <<- TRUE
                    datCounter <<- datCounter + 1
                }
            }, text = function(text) {
                if (datBool) {
                    parsedXML[[rowCounter]][[colCounter]][datCounter] <<- text
                }
            }, endElement = function(name, uri) {
                if (name == "COL") {
                    colBool <<- FALSE
                    datCounter <<- 0
                }
                if (name == "DATA") {
                    datBool <<- FALSE
                }
            }
        ))
    return(do.call(rbind, parsedXML))
    }

########################################################################################
# Input: File 1 xml path, temporary folder where raw files will be copied and converted,
#  path to msconvert.exe from ProteoWizard
# Output: Parsed file 1 (list), Parsed file 1 (matrix), minimum and maximum column 
#  numbers, minimum and maximum fraction numbers
# 
# This method parses file 1 xml into two data structures, one is a list and another is
#  a matrix. Both are used later on in findEmptyPeaks. First the xml is transformed into
#  a list with the XML package, so empty tag values are retained. Because each tag that
#  holds data has the same name, indices had to be used to access their data. Column 
#  numbers, fraction numbers, and paths are retained, along with a xcmsRaw filler space
#  that will be later populated with the xcmsRaw object associated with the path at that
#  row. Then each row is then named the path that is in that row. Minimum and maximum 
#  column and fractions are calculated, which are later used in findEmptyPeaks as a 
#  range where data will be looked at.
# A matrix is then also created, with the path of the object at the column number (row)
#  and fraction number (column) as the indices. Both a list and a matrix are created so
#  searching for the xcmsRaw object by indices and name (path) would be possible instead
#  of looking through the matrix for the path to get the column and fraction numbers 
#  (old design). There is also a check for if there are multiple paths at the same 
#  column and fraction number pair.
# With the list, for each row, the raw file at the path provided is then copied over to
#  the temporary folder, and converted to an mzXML file in the same directory. Then the
#  mzXML file is fed into xcms to be converted to an xcmsRaw object, and stored into the
#  list.
# https://xcmsonline.scripps.edu/docs/fileformats.html
########################################################################################
parseFile1 <- function(file1Path, columnLength) {
    print(paste("Parsing ", file1Path, sep=""))
    # if (!file.exists(msconv)) waitForInput(paste(msconv," does not exist. Please check ProteoWizard installation and config file.", sep=""))
    if (!file.exists(file1Path)) waitForInput(paste(file1Path, " does not exist. Please check input folder.", sep=""))

    parsedXML <- parseFMXML(file1Path)
    parsedF1 <- parsedXML
    parsedF1[,1:2] <- as.numeric(parsedF1[,1:2])
    colnames(parsedF1) <- c("column", "fraction", "path")
    rownames(parsedF1) <- c(unlist(parsedF1[,"path"]))

    columns <- unlist(parsedF1[,"column"])
    fractions <- unlist(parsedF1[,"fraction"])

    parsedF1 <- parsedF1[order(columns, fractions),]

    minCol <- min(columns)
    maxCol <- max(columns)
    minFrac <- min(fractions)
    maxFrac <- max(fractions)

    return(list(parsedF1=parsedF1, minCol=minCol, maxCol=maxCol, minFrac=minFrac, maxFrac=maxFrac))
}

########################################################################################
# Input: File 2 xml path, SILAC boolean
# Output: Parsed file 2
# 
# This method parses file 2 xml into a list. Depending if the data is SILAC or
#  LabelFree, file 2 xml is different, and is parsed differently. For example, the 
#  ordering of the data SILAC and LabelFree differ in where the raw file paths are, and 
#  LabelFree has an extra field for standard deviation. Also, the number of masses (2 
#  for SILAC, 1 for LabelFree), will affect the output and what is being stored. 
#  Essentially there are peak areas, peak widths, and metric (scores) for every mass and
#  retention time (aligned) combination. Peptide counter, retention times (vector), 
#  paths (vector), mass (vector), average retention time (number), aligned retention 
#  times (vector), fraction number (vector), peak areas (list(vector) *# mass), peak 
#  widths (list(vector) *# mass), and metric scores (list(vector) *# mass) are stored
#  per row in the output list, with information associated with file 2 xml.
# Changes: Now also parses actual RT's from file 2 (new xml format, extra repeating
#  field. comma-delimited RTs, one per replicate).
########################################################################################
parseFile2 <- function(file2Path, SILAC) {
    print(paste("Parsing ", file2Path, sep=""))
    if (!file.exists(file2Path)) waitForInput(paste(file2Path, " does not exist. Please check input folder.", sep=""))

    parsedXML <- parseFMXML(file2Path)
    counterIndex <- 1
    rtIndex <- 2
    pathIndex <- ifelse(SILAC, 3, 5)
    massIndex <- 4
    if (SILAC) {
        actRtIndex <- c(5,6,7)
    } else {
        actRtIndex <- 6
    }
    repSize <- max(unlist(lapply(parsedXML[1, actRtIndex[1]][[1]], function(str) length(gregexpr(",", str)[[1]]) + 1)))

    parsedXML[, counterIndex] <- as.numeric(parsedXML[, counterIndex])
    parsedXML[, rtIndex] <- lapply(1:nrow(parsedXML), function(j) {
        rtL <- as.numeric(unlist(parsedXML[j, rtIndex]))
        rtL[is.na(rtL)] <- 0
        return(rtL)
        })

    # This will, for each peptide/row, parse the actual retention times 
    #  for all replicates and columns provided. The outer lapply goes through
    #  all peptides/rows, and the inner lapply goes through all columns. For
    #  each column it will split the string of retention times by "," and 
    #  fill in blanks with "NA", into a list. rbind is called to flatten 
    #  all the lists and merge them into a matrix.
    for (num in actRtIndex) {
        parsedXML[, num] <- lapply(1:nrow(parsedXML), function(j) {
            l <- lapply(1:length(parsedXML[j, num][[1]]), function(k) {
                x <- as.numeric(unlist(strsplit(unlist(as.character(parsedXML[j, num][[1]][k])), ",")))
                length(x) <- repSize
                x[is.na(x)] <- 0
                return(x)
                })
            do.call(rbind, l)
            }) 
    }

    if (SILAC) {
        parsedXML[, massIndex] <- lapply(1:nrow(parsedXML), function(i) 
            as.numeric(unlist(strsplit(unlist(parsedXML[i, massIndex]), ","))))
        colnames(parsedXML) <- c("counter", "rt", "path", "mass", "actRt1", "actRt2", "actRt3")
    } else {
        parsedXML[, massIndex] <- as.numeric(parsedXML[, massIndex])
        parsedXML <- parsedXML[,c(counterIndex, rtIndex, pathIndex, massIndex, actRtIndex)]
        colnames(parsedXML) <- c("counter", "rt", "path", "mass", "actRt")
    }

    columnLength <- length(unlist(parsedXML[1,rtIndex]))
    massLength <- length(unlist(parsedXML[1,massIndex]))
    # avRT <- lapply(1:nrow(parsedXML), function(k) mean(unlist(parsedXML[k, rtIndex]), na.rm=TRUE))
    # alRT <- list(as.numeric(rep(0, columnLength)))
    # fraction <- if(SILAC) list(as.numeric(rep(-1, columnLength))) else list(as.numeric(rep(0, columnLength)))
    # peakArea <- list(rep(list(rep("", columnLength)), each=massLength))
    # peakWidth <- list(rep(list(rep(0, columnLength)), each=massLength))
    # metric <- list(rep(list(rep("no score", columnLength)), each=massLength))

    # parsedF2 <- cbind(parsedXML, avRT, alRT, fraction, peakArea, peakWidth, metric)

    parsedF2 <- parsedXML
    npeptides <- nrow(parsedF2)

    return(list(parsedF2=parsedF2, massLength=massLength, columnLength=columnLength, npeptides=npeptides))
}

########################################################################################
# Input: align.txt file path
# Output: Parsed align file
# 
# This method parses the align text file from quickPeakAlignment as a table, separated
#  by tabs. The first column string (3 numbers separated by ".") is stored under "a", 
#  and the second column string (same format as first column) is stored under "b". Then 
#  the fifth column of retention times tuples are split by ";", and the first value is 
#  stored under "RTa" and the second value is stored under "RTb", both as vectors. This 
#  is done for every row in the align text file.
# Changes: Now filters rt so that both RTa and RTb are increasing (smoothing out data).
#  This is done by finding "peaks" and "valleys", and the edges are handled separately. 
#  Comparison between each list of numbers and its sorted form are made to check if they 
#  are in order. The way peaks and valleys in the retention time vectors are found is by
#  taking the differences of the numbers and finding the sign (-1 and 1 for negative
#  and positive numbers respectively) of them. Then the difference of them is taken
#  again, and if there were any non-zero numbers (which will be only -2 and 2), then
#  they are removed. Inf and -Inf are padded to the ends of the vectors to consider the
#  edges of the vectors. This will be added to quickPeakAlignment.
########################################################################################
parseAlign <- function(alignPath) {
    print(paste("Parsing ", alignPath, sep=""))
    if (!file.exists(alignPath)) waitForInput(paste(alignPath, " does not exist. Please check input folder.", sep=""))
    alignFile <- read.table(alignPath, sep="\t", fileEncoding="UTF-16")
    nRowAlign <- nrow(alignFile)
    parsedAlign <- matrix(list(), nrow=nRowAlign, ncol=4)
    for (i in 1:nRowAlign) {
        parsedAlign[i, 1] <- list(as.character(gsub("\\[", "", alignFile[i,][["V1"]])))
        parsedAlign[i, 2] <- list(as.character(alignFile[i,][["V2"]]))

        v5 <- unlist(strsplit(gsub("\\]", "", as.character(alignFile[i,][["V5"]])), ";"))
        v5a <- as.numeric(lapply(v5, function(a) as.numeric(gsub("\\(", "", unlist(strsplit(a, ","))[1]))))
        v5b <- as.numeric(lapply(v5, function(b) as.numeric(gsub("\\)", "", unlist(strsplit(b, ","))[2]))))

        if (v5a != sort(v5a) || v5b != sort(v5b)) {
            matchv5av <- which(diff(sign(diff(c(Inf, v5a, Inf))))==2)
            if (1 %in% matchv5av) {
                if (v5a[1] <= v5a[2]) {
                    matchv5av <- matchv5av[matchv5av != 1]
                }
            }
            v5av <- v5a[-matchv5av]
            if (length(v5av) == 0) v5av <- v5a
            matchv5ap <- which(diff(sign(diff(c(-Inf, v5av, -Inf))))==-2)
            if (length(v5av) %in% matchv5ap) {
                if (v5av[length(v5av)] >= v5av[length(v5av)-1]) {
                    matchv5ap <- matchv5ap[matchv5ap != length(v5av)]
                }
            }
            v5ap <- v5av[-matchv5ap]
            if (length(v5ap) == 0) v5ap <- v5av

            if (length(v5ap) == 0) {
                if (length(v5av) == 0) {
                    v5atmp <- v5a
                } else {
                    v5atmp <- v5av
                }
            } else {
                v5atmp <- v5ap
            }

            matchv5bv <- which(diff(sign(diff(c(Inf, v5b, Inf))))==2)
            if (1 %in% matchv5bv) {
                if (v5b[1] <= v5b[2]) {
                    matchv5bv <- matchv5bv[matchv5bv != 1]
                }
            }
            v5bv <- v5b[-matchv5bv]
            if (length(v5bv) == 0) v5bv <- v5b
            matchv5bp <- which(diff(sign(diff(c(-Inf, v5bv, -Inf))))==-2)
            if (length(v5bv) %in% matchv5bp) {
                if (v5bv[length(v5bv)] >= v5bv[length(v5bv)-1]) {
                    matchv5bp <- matchv5bp[matchv5bp != length(v5bv)]
                }
            }
            v5bp <- v5bv[-matchv5bp]
            if (length(v5bp) == 0) v5bp <- v5bv

            if (length(v5bp) == 0) {
                if (length(v5bv) == 0) {
                    v5btmp <- v5b
                } else {
                    v5btmp <- v5bv
                }
            } else {
                v5btmp <- v5bp
            }

            matcha <- match(v5atmp, v5a)
            matchb <- match(v5btmp, v5b)
            match <- intersect(matcha[!is.na(matcha)], matchb[!is.na(matchb)])

            v5a <- v5a[match]
            v5b <- v5b[match]
        }

        parsedAlign[i, 3] <- list(v5a)
        parsedAlign[i, 4] <- list(v5b)
    }

    colnames(parsedAlign) <- c("a", "b", "RTa", "RTb")
    return(parsedAlign)
}

########################################################################################
# Input: xcal_RT_from_rep.xml file path
# Output: Parsed xcal xml file
# 
# This method parses the xcal xml as a list, to retain empty tag values in the XML 
#  (NA). For every row, the peptide counter and replicate number (vector) of each column
#  are stored.
########################################################################################
parsexCal <- function(xCalPath) {
    print(paste("Parsing ", xCalPath, sep=""))
    if (!file.exists(xCalPath)) waitForInput(paste(xCalPath, " does not exist. Please check input folder.", sep=""))
    parsedXML <- parseFMXML(xCalPath)
    colnames(parsedXML) <- c("counter", "rep")
    parsedXML[, "counter"] <- as.numeric(parsedXML[, "counter"])
    parsedXML[, "rep"] <- lapply(1:nrow(parsedXML), function(j) as.numeric(unlist(parsedXML[j, "rep"])))

    return(parsedXML)
}

########################################################################################
# Input: ta (key/string of "column.fraction.replicate number from xcal xml"), tb 
#  (key/string of "column.fraction.replicate number of input file xml"), retention time
#  that needs to be adjusted, parsedAlign list
# Output: Aligned retention time (adjusted)
# 
# This method takes in a retention time and returns a aligned retention time, using data
#  from quickPeakAlignment. It goes through the parsed list from quickPeakAlignment 
#  until it finds a match in both ta and tb. Then depending on the specific case, 
#  indices are specified that are later used in the calculation of the aligned retention
#  time. If the retention time is less than the first value in RTa, then the first two 
#  numbers in both RTa and RTb are used. If the retention time is greater than the last 
#  value in RTa, then the last two numbers in both RTa and RTb are used. If neither 
#  cases are met, then it goes through the list of values in RTa until the retention 
#  time provided is between two values in RTa.
# Finally, the aligned retention time is calculated using Euler's method of 
#  approximation.
# An edge case is considered when there are no matches, or if the indices where there 
#  is a match are the same, then the original retention time provided will be returned
#  instead.
# This method is a copy of the method from the original AutoFill (visual basic)
########################################################################################
alignRT <- function(ta, tb, rt, parsedAlign) {
    match <- intersect(which(parsedAlign[,"a"] == ta), which(parsedAlign[,"b"] == tb))
    if (length(match) > 0 && min(match) != Inf) {
        i <- min(match)
        p1 <- 1
        p2 <- 1
        j1 <- 1
        j2 <- length(parsedAlign[i,"RTa"][[1]])
        if (rt < parsedAlign[i,"RTa"][[1]][j1]) {
            p1 <- 1
            p2 <- 2
        } else if (rt >= parsedAlign[i,"RTa"][[1]][j2]) {
            p1 <- j2 - 1
            p2 <- j2
        } else {
            for (j in j1:(j2-1)) {
                if (rt >= parsedAlign[i,"RTa"][[1]][j] && rt < parsedAlign[i,"RTa"][[1]][j + 1]) {
                    p1 <- j
                    p2 <- j + 1
                }
            }
        }
        if (p1 == p2) {
            return(rt)
        } else {
            return((parsedAlign[i,"RTb"][[1]][p2] - parsedAlign[i,"RTb"][[1]][p1]) 
                / (parsedAlign[i,"RTa"][[1]][p2] - parsedAlign[i,"RTa"][[1]][p1]) 
                * (rt - parsedAlign[i,"RTa"][[1]][p1]) + parsedAlign[i,"RTb"][[1]][p1])
        }
    }
    return(rt)
}

########################################################################################
# Input: Peptide counter (tcounter), column number, parsed xcal
# Output: Replicate number associated with peptide and column
# 
# This method goes through the parsed xcal until it finds the peptide. From there, it 
#  returns the replicate number associated with that peptide and the column number that
#  was provided.
########################################################################################
getRepNumberOfThisRT <- function(tcounter, col, xCal) {
    match <- which(xCal[,"counter"] == tcounter)
    if (length(match) > 0) {
        return(xCal[min(match),"rep"][[1]][col])
    }
    return("")
}

########################################################################################
# Input: Row of parsed file 2, parsed file 1 list, parsed file 1 matrix, minimum column
#  number, maximum column number, minimum fraction number, maximum fraction number, 
#  mass window (from config file, aka ppm), replicate number (based on file name and 
#  supplied by run.r), parsed align file, parsed xcal file, svm object for metrics, 
#  SILAC boolean
# Output: Filled row
# 
# This method is the central part of AutoFill. First you take a row and find the 
#  average of the retention times. If the retention time average is less than 0, then 
#  skip the row. This can occur when there are no retention times supplied or if all of
#  the retention times are zero. Skip the row also if the mass supplied are less than 0.
# Now this will iterate through the columns of the row, and calculates an aligned 
#  retention time with the results from quickPeakAlignment and information from the run.
#  If there is a retention time and raw file path for the column of the peptide, then it
#  will do the following:
#   The first key is based on the column number, fraction number, and replicate number
#   provided in xcal xml for the associated peptide and column number. Key 2 is the same
#   but with a different replicate number, the one from the input (based on the input 
#   files). They are then used to generate an aligned retention time using alignRT. 
#   Then with the aligned retention time, peakCalc will run from xcms2.r, using xcms to 
#   find a peak (details can be found there). If there is a peak returned, then it will 
#   be added to the output of the row.
#  If something is missing in the current column, it will try to do the following:
#   Using the other columns of the same peptide, it will find the aligned retention time
#   in the same manner as before, but use the retention time of the other columns, and 
#   run peakCalc with it. Now there are a few cases it needs to consider. Of the results
#   with the information from the other columns, it will try to retain the peak area 
#   and associated information from a detected peak, and not a calculated peak. If no 
#   peaks were detected but instead calculated, it will take the one with the largest
#   peak area (intensity). This is also true for the case where multiple peaks were 
#   detected.
# Changes from the original: No fraction support (original had that part commented out),
#  gives priority to detected peaks vs adding area in region.
########################################################################################
findEmptyPeaks <- function(row, column, xcmsRaw, peak, parsedF1, minCol, maxCol, minFrac, maxFrac, massWindow, repNum, parsedAlign, 
    parsedxCal, svmObject=NULL, SILAC) {

    #print(row[["counter"]])

    if (is.na(row[["mass"]]) || any(row[["mass"]] < 0)) {
        return(peak)
    }
    if (SILAC) {
        if (row[["rt"]][column] > 0 && !is.na(row[["path"]][column])) {
            myfrac <- parsedF1[row[["path"]][column],"fraction"][[1]]
            key1 <- paste(column, myfrac, getRepNumberOfThisRT(row[["counter"]], column, parsedxCal), sep=".")
            key2 <- paste(column, myfrac, repNum, sep=".")
            tmpRT <- alignRT(key1, key2, row[["rt"]][column], parsedAlign)
            p <- peakCalc(xcmsRaw, row[["mass"]], tmpRT*60, massWindow, svmObject, TRUE, FALSE)
            for (l in 1:length(p[["peaks"]])) {
                if(!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                    peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                    peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                    peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                    peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                    peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                    if (p[["peaks"]][[l]][["metric"]][["metric"]] != "no score") {
                        peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                        peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                        peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                        peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
						peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
						peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
						peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
						peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
						peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
						peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
						peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
						peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
						peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                    }
                }
            }
            peak[,"fraction"] <- myfrac
            peak[,"alRT"] <- tmpRT
        } else {
            maxArea <- -Inf
            for (k in minCol:maxCol) {
                if (row[["rt"]][k] > 0 && !is.na(row[["path"]][k])) {
                    myfrac <- parsedF1[row[["path"]][k],"fraction"][[1]]
                    key1 <- paste(k, myfrac, getRepNumberOfThisRT(row[["counter"]], k, parsedxCal), sep=".")
                    key2 <- paste(column, myfrac, repNum, sep=".")
                    tmpRT <- alignRT(key1, key2, row[["rt"]][k], parsedAlign)
                    p <- peakCalc(xcmsRaw, row[["mass"]], tmpRT*60, massWindow, svmObject, TRUE, FALSE)
                    peakAreas <- unlist(lapply(1:length(p[["peaks"]]), function(n) p[["peaks"]][[n]][["peakArea"]]))
                    # Edge case if peak areas returned contain "NA" or are only "NA"
                    # Min is used if all peak areas returned are negative
                    suppressWarnings({
                        if (any(peakAreas > 0, na.rm=TRUE)) {
                            tmpMaxArea <- max(peakAreas, na.rm=TRUE)
                        } else {
                            tmpMaxArea <- min(peakAreas, na.rm=TRUE)
                        }
                        })
                    # If max of a vector of "NA"'s or min of a vector of "NA"'s is called
                    #  They will return -Inf or Inf, respectively
                    if (tmpMaxArea != -Inf && tmpMaxArea != Inf) {
                        if (maxArea == -Inf) { # Initial case, put in a peak if found and is first one.
                            maxArea <- tmpMaxArea
                            for (l in 1:length(p[["peaks"]])) {
                                if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                    peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                    peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                    peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                    peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                    peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                    peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                    peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                    peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                    peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
									peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
									peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
									peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
									peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
									peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
									peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
									peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
									peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
									peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                                }
                            }
                            peak[,"fraction"] <- myfrac
                            peak[,"alRT"] <- tmpRT
                        } else {
                            if (tmpMaxArea * maxArea > 0) { # Same sign comparision. Instead of having a if/else if, combining both by
                                                            #  multiplying them and checking that sign will tell whether the peak areas
                                                            #  were of the same sign. Then the comparision of which is greater (more 
                                                            #  negative or positive, depending on the sign of the ares) is made.
                                if (abs(tmpMaxArea) > abs(maxArea)) {
                                    maxArea <- tmpMaxArea
                                    for (l in 1:length(p[["peaks"]])) {
                                        if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                            peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                            peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                            peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                            peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                            peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                            peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                            peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                            peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                            peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
											peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
											peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
											peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
											peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
											peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
											peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
											peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
											peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
											peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                                        }
                                    }
                                    peak[,"fraction"] <- myfrac
                                    peak[,"alRT"] <- tmpRT
                                }
                            } else {
                                if (tmpMaxArea > 0 && maxArea < 0) { # Give priority to detected peak (positive) over calculated peak (negative)
                                    maxArea <- tmpMaxArea
                                    for (l in 1:length(p[["peaks"]])) {
                                        if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                            peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                            peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                            peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                            peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                            peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                            peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                            peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                            peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                            peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
											peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
											peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
											peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
											peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
											peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
											peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
											peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
											peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
											peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                                        }
                                    }
                                    peak[,"fraction"] <- myfrac
                                    peak[,"alRT"] <- tmpRT
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (row[["actRt"]][column,rep] > 0 && !is.na(row[["path"]][column])) {
            myfrac <- parsedF1[row[["path"]][column],"fraction"][[1]]
            tmpRT <- row[["actRt"]][column,rep]
            p <- peakCalc(xcmsRaw, row[["mass"]], tmpRT*60, massWindow, svmObject, TRUE, FALSE)
			for (l in 1:length(p[["peaks"]])) {
                if(!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                    peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                    peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                    peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                    peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                    peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                    peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                    peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                    peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                    peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
					peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
					peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
					peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
					peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
					peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
					peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
					peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
					peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
					peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
					peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                }
            }
            peak[,"fraction"] <- myfrac
            peak[,"alRT"] <- tmpRT
        } else {
            maxArea <- -Inf
            for (k in minCol:maxCol) {
                for (m in 1:length(row[["actRt"]][k,])) {
                    if (row[["actRt"]][k,m] > 0 && !is.na(row[["path"]][k])) {
                        myfrac <- parsedF1[row[["path"]][k],"fraction"][[1]]
                        key1 <- paste(k, myfrac, m, sep=".")
                        key2 <- paste(column, myfrac, repNum, sep=".")
                        tmpRT <- alignRT(key1, key2, row[["actRt"]][k,m], parsedAlign)
                        p <- peakCalc(xcmsRaw, row[["mass"]], tmpRT*60, massWindow, svmObject, TRUE, FALSE)
                        peakAreas <- unlist(lapply(1:length(p[["peaks"]]), function(n) p[["peaks"]][[n]][["peakArea"]]))
                        # Edge case if peak areas returned contain "NA" or are only "NA"
                        # Min is used if all peak areas returned are negative
                        suppressWarnings({
                            if (any(peakAreas > 0, na.rm=TRUE)) {
                                tmpMaxArea <- max(peakAreas, na.rm=TRUE)
                            } else {
                                tmpMaxArea <- min(peakAreas, na.rm=TRUE)
                            }
                            })
                        # If max of a vector of "NA"'s or min of a vector of "NA"'s is called
                        #  They will return -Inf or Inf, respectively
                        if (tmpMaxArea != -Inf && tmpMaxArea != Inf) {
                            if (maxArea == -Inf) { # Initial case, put in a peak if found and is first one.
                                maxArea <- tmpMaxArea
                                for (l in 1:length(p[["peaks"]])) {
                                    if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                        peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                        peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                        peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                        peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                        peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                        peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                        peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                        peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                        peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
										peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
										peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
										peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
										peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
										peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
										peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
										peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
										peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
										peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                                    }
                                }
                                peak[,"fraction"] <- myfrac
                                peak[,"alRT"] <- tmpRT
                            } else {
                                if (tmpMaxArea * maxArea > 0) { # Same sign comparision. Instead of having a if/else if, combining both by
                                                                #  multiplying them and checking that sign will tell whether the peak areas
                                                                #  were of the same sign. Then the comparision of which is greater (more 
                                                                #  negative or positive, depending on the sign of the ares) is made.
                                    if (abs(tmpMaxArea) > abs(maxArea)) {
                                        maxArea <- tmpMaxArea
                                        for (l in 1:length(p[["peaks"]])) {
                                            if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                                peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                                peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                                peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                                peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                                peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                                peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                                peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                                peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                                peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
												peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
												peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
												peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
												peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
												peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
												peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
												peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
												peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
												peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
                                            }
                                        }
                                        peak[,"fraction"] <- myfrac
                                        peak[,"alRT"] <- tmpRT
                                    }
                                } else {
                                    if (tmpMaxArea > 0 && maxArea < 0) { # Give priority to detected peak (positive) over calculated peak (negative)
                                        maxArea <- tmpMaxArea
                                        for (l in 1:length(p[["peaks"]])) {
                                            if (!is.na(p[["peaks"]][[l]][["peakArea"]])) {
                                                peak[l,"peakRTL"] <- p[["peaks"]][[l]][["peakLeftRT"]]
                                                peak[l,"peakRT"] <- p[["peaks"]][[l]][["peakApexRT"]]
                                                peak[l,"peakRTR"] <- p[["peaks"]][[l]][["peakRightRT"]]
                                                peak[l,"peakArea"] <- p[["peaks"]][[l]][["peakArea"]]
                                                peak[l,"peakWidth"] <- p[["peaks"]][[l]][["peakWidth"]]
                                                peak[l,"metric"] <- p[["peaks"]][[l]][["metric"]][["metric"]]
                                                peak[l,"sn"] <- p[["peaks"]][[l]][["sn"]]
                                                peak[l,"kurtosis"] <- p[["peaks"]][[l]][["metric"]][["kurtosis"]]
                                                peak[l,"skew"] <- p[["peaks"]][[l]][["metric"]][["skew"]]
												peak[l,"derivativeCount"] <- p[["peaks"]][[l]][["metric"]][["derivativeCount"]]
												peak[l,"ms2Call"] <- p[["peaks"]][[l]][["metric"]][["ms2Call"]]
												peak[l,"zeroCount"] <- p[["peaks"]][[l]][["metric"]][["zeroCount"]]
												peak[l,"maxIntensity"] <- p[["peaks"]][[l]][["metric"]][["maxIntensity"]]
												peak[l,"rtWindow"] <- p[["peaks"]][[l]][["metric"]][["rtWindow"]]
												peak[l,"peakSig"] <- p[["peaks"]][[l]][["metric"]][["peakSig"]]
												peak[l,"peakSharp"] <- p[["peaks"]][[l]][["metric"]][["peakSharp"]]
												peak[l,"TPASR"] <- p[["peaks"]][[l]][["metric"]][["TPASR"]]
												peak[l,"FWHM"] <- p[["peaks"]][[l]][["metric"]][["FWHM"]]
											}
                                        }
                                        peak[,"fraction"] <- myfrac
                                        peak[,"alRT"] <- tmpRT
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    # print(peak)
    return(list(peak))
}

########################################################################################
# Input: Output from findEmptyPeaks (merged), output file path, SILAC boolean
# Output: Nothing
# 
# This method writes the output from findEmptyPeaks into a form that filemaker will be
#  parsing. The formats are based on the visual basic AutoFill output format. Depending
#  on if the data is SILAC or LabelFree, the output format differs. 
# SILAC format: 
#  Peptide counter    Aligned RT1|Fraction: peak area, peakWidth, score (of mass 1)|
#                                           peak area, peakWidth, score (of mass 2)
#                     Aligned RT2|Fraction: ... 
# LabelFree format:
#  Peptide counter    Aligned RT1, peak area, peak width, fraction, score    Aligned RT2
#                     ...
########################################################################################
writeOutput <- function(parsedF2, resultsF2, massLength, columnLength, npeptides, outputPath, SILAC) {
    print("Buffering output.")
    f <- file(outputPath)
    counters <- unlist(parsedF2[,"counter"])
    if (SILAC) {
        # lines <- lapply(resultsF2, function(row) {
        #     output <- paste(as.character(row[["counter"]]), "\t", sep="")
        #     for (j in 1:length(row[["rt"]])) {
        #         output <- paste(output, row[["alRT"]][j], "|", row[["fraction"]][j], ": ", sep="")
        #         for (k in 1:length(row[["mass"]])) {
        #             output <- paste(output, row[["peakArea"]][[k]][j], ",", row[["peakWidth"]][[k]][j], ",", 
        #                 row[["metric"]][[k]][j], "|", sep="")
        #             if (k == (length(row[["mass"]]))) {
        #                 output <- paste(output, "\t", sep="")
        #             }
        #         }
        #     }
        #     return(output)
        #     })
        lines <- lapply(1:npeptides, function(i) {
            output <- paste(counters[i], "\t", sep="")
            for (j in 1:columnLength) {
                output <- paste(output, resultsF2[[i,j]][1,"alRT"], "|", resultsF2[[i,j]][1,"fraction"], ": ", sep="")
                for (k in 1:massLength) {
#                    output <- paste(output,",",resultsF2[[i,j]][1,"peakRTL"],",",resultsF2[[i,j]][1,"peakRT"],",",resultsF2[[i,j]][1,"peakRTR"], resultsF2[[i,j]][k,"peakArea"], ",", resultsF2[[i,j]][k,"peakWidth"], ",", 
#                        resultsF2[[i,j]][k,"metric"],",",resultsF2[[i,j]][k,"sn"],",",resultsF2[[i,j]][k,"kurtosis"],",",
#                        resultsF2[[i,j]][k,"skew"],"|", sep="")
                    output <- paste(output, resultsF2[[i,j]][k,"peakArea"], ",",
					  resultsF2[[i,j]][k,"peakWidth"],",", 
                      resultsF2[[i,j]][k,"metric"],",",
					  resultsF2[[i,j]][k,"sn"],",",
					  resultsF2[[i,j]][k,"kurtosis"],",",
                      resultsF2[[i,j]][k,"skew"],",",	#everything below this line added 10242016
					  resultsF2[[i,j]][k,"peakRT"],",",
					  resultsF2[[i,j]][k,"derivativeCount"],",",
					  resultsF2[[i,j]][k,"ms2Call"],",",
					  resultsF2[[i,j]][k,"zeroCount"],",",
					  resultsF2[[i,j]][k,"maxIntensity"],",",
					  resultsF2[[i,j]][k,"rtWindow"],",",
					  resultsF2[[i,j]][k,"peakSig"],",",
					  resultsF2[[i,j]][k,"peakSharp"],",",
					  resultsF2[[i,j]][k,"TPASR"],",",
					  resultsF2[[i,j]][k,"FWHM"],
					  "|", sep="")
                }
                if (j < columnLength) {
                    output <- paste(output, "\t", sep="")
                }
            }
            return(output)
            })
    } else {
        lines <- lapply(1:npeptides, function(i) {
            output <- paste(counters[i], "\t", sep="")
            for (j in 1:columnLength) {
                # output <- paste(output, resultsF2[[i,j]][1,"alRT"], ",",resultsF2[[i,j]][1,"peakRTL"], ",",
				    # resultsF2[[i,j]][1,"peakRT"], ",",resultsF2[[i,j]][1,"peakRTR"], ",", resultsF2[[i,j]][1,"peakArea"], ",", 
                    # resultsF2[[i,j]][1,"peakWidth"], ",", resultsF2[[i,j]][1,"fraction"], ",", 
                    # resultsF2[[i,j]][1,"metric"], ",",resultsF2[[i,j]][1,"sn"], ",",resultsF2[[i,j]][1,"kurtosis"], ",",
                    # resultsF2[[i,j]][1,"skew"],sep="")
                output <- paste(output, resultsF2[[i,j]][1,"alRT"],",",
				  resultsF2[[i,j]][1,"peakArea"],",", 
                  resultsF2[[i,j]][1,"peakWidth"],",",
				  resultsF2[[i,j]][1,"fraction"],",", 
                  resultsF2[[i,j]][1,"metric"],",",
				  resultsF2[[i,j]][1,"sn"],",",
				  resultsF2[[i,j]][1,"kurtosis"],",",
                  resultsF2[[i,j]][1,"skew"],",",
				  resultsF2[[i,j]][1,"peakRT"],",", #everything below this line added 10242016
				  resultsF2[[i,j]][1,"derivativeCount"],",",
				  resultsF2[[i,j]][1,"ms2Call"],",",
				  resultsF2[[i,j]][1,"zeroCount"],",",
				  resultsF2[[i,j]][1,"maxIntensity"],",",
				  resultsF2[[i,j]][1,"rtWindow"],",",
				  resultsF2[[i,j]][1,"peakSig"],",",
				  resultsF2[[i,j]][1,"peakSharp"],",",
				  resultsF2[[i,j]][1,"TPASR"],",",
				  resultsF2[[i,j]][1,"FWHM"],
				  sep="")
                if (j < columnLength) {
                    output <- paste(output, "\t", sep="")
                }
            }
            return(output)
            })
    }
    print(paste("Writing output to: ", outputPath, sep=""))
    writeLines(unlist(lines), f)
    close(f)
}

########################################################################################
# Input: Config file path
# Output: Parsed config file list
# 
# This method parses the config file generated by filemaker and splits by ";". The first
#  value is the key (name), and the second value is the value of the key. This 
#  is basically a dictionary in R. File names, constants, etc. are in this config file 
#  and parsed config file list.
########################################################################################
constants <- function(f) {
    print(paste("Parsing ", f, sep=""))
    parsedConstants <- list()
    cfgFile <- file(f, open="r", encoding="UTF-16LE")
    lines <- readLines(cfgFile, warn=FALSE)
    close(cfgFile)
    for (line in lines) {
        split <- unlist(strsplit(line, ";"))
        parsedConstants[[split[1]]] <- split[2]
    }
    return(parsedConstants)
}

########################################################################################
# Input: String (message)
# Output: Nothing
# 
# This method prompts the user when it is called, providing the option to continue or 
#  quit the current running script, based on [key press] and [enter]. It also will try
#  again if an improper key was entered.
# http://stackoverflow.com/questions/10266963/moving-files-between-folders
########################################################################################
waitForInput <- function(error="") {
    cat(paste("Error occured: ", error, "\n", sep=""))
    cat("To continue, type [y]. To cancel, type [n]. Then press [enter].\n")
    key <- scan("stdin", n=1, what="")
    if (key == "y") {
        cat("Continuing job.\n")
    } else if (key == "n") {
        cat("Canceling job.\n")
        quit()
    } else {
        cat("Invalid key input, try again\n")
        waitForInput()
    }
}

########################################################################################
# Input: Nothing
# Output: Directory of running script
# 
# This method looks at the argument, tries to get the path of the current running 
#  script, and return it.
# http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
########################################################################################
currentDir <- function() {
    initial.options <- commandArgs(trailingOnly = FALSE)
    file.arg.name <- "--file="
    script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
    script.basename <- gsub("\\\\", "/", dirname(script.name))
    return(script.basename)
}

########################################################################################
# Input: Directory path
# Output: Nothing
# 
# This method will remove all files in the directory provided.
########################################################################################
removeOldFiles <- function(dir) {
    oldFiles <- list.files(dir)
    unlink(lapply(oldFiles, function(names) paste(dir, names, sep="/")), recursive = TRUE, force = FALSE)
}


########################################################################################


########################################################################################
# Workflow:
#  Get directory of running script
#  Source external scripts in same directory
#  Parse arguments:
#   1. Config file path
#   2. Replicate number
#   3. SILAC or LABELFREE indicator
#   4. PARALLEL processing indicator
#  Form paths of input data
#  Parse File 1 xml
#  Parse File 2 xml
#  Parse align.txt
#  Parse xcal_RT_from_rep.xml
#  Parse trainingSet.csv and create svm for quality peak metrics
#  Detect number of cores, and run findEmptyPeaks for every row in file 2/peptide, 
#   either in one process or multiple processesm and update progress bar
#  Write output to file
# 
# Progress bar design: There is a temporary file that for each run of findEmptyPeaks 
#  (for each row in file 2/peptide), it will write a "0", and count the number of them.
#  Based on the number of them, the progress bar will be updated.
# 
# Multiprocessing design: This uses the built in parallel package bundled with R 2.14+.
#  The current design uses parts of SNOW, and creates a cluster. The parsed objects are 
#  exported to each of the processes in the cluster. Then findEmptyPeaks will run in the 
#  same manner as with one process, but with a parallel version of apply. Also file 1 
#  from the master process is removed before quantitation and the garbage collector is 
#  invoked to remove it, freeing up memory early on once the parsed file 1 is copied 
#  over to the processes in the cluster. A progress bar is implemented in that it will
#  create blank files based on the peptide counter and count them.
########################################################################################
ptm <- proc.time()

scriptDir <- currentDir()
if (length(scriptDir) == 0) {
    scriptDir <- "."
}
xcmsScriptPath <- paste(sep="/", scriptDir, xcmsScript)
cat(paste("Sourcing", xcmsScriptPath, "from", xcmsScriptPath, "\n", sep=" "))
source(xcmsScriptPath)

args <- commandArgs(trailingOnly = TRUE)
config <- args[1]
rep <- as.numeric(args[2])

parsedConstants <- constants(config)

if (toupper(args[3]) == "SILAC") {
    SILAC <- TRUE
    file1 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], rep, parsedConstants[["f1s"]], sep="")
    file2 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], rep, parsedConstants[["f2s"]], sep="")
    outputPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["coms"]], rep, ".txt", sep="")
    skippedPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["skips"]], rep, ".txt", sep="")
} else if (toupper(args[3]) == "LABELFREE") {
    SILAC <- FALSE
    file1 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], rep, parsedConstants[["f1lf"]], sep="")
    file2 <- paste(parsedConstants[["xmldir"]], parsedConstants[["rep"]], rep, parsedConstants[["f2lf"]], sep="")
    outputPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["comlf"]], rep, ".txt", sep="")
    skippedPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["skiplf"]], rep, ".txt", sep="")
} else {
    stop("Please provide type of data for quantitation.")
}

if (toupper(args[4]) == "PARALLEL") {
    cores <- detectCores()
} else if (toupper(args[4]) == "NOPARALLEL") {
    cores <- 1
} else {
    stop("Please provide multiprocessing option")
}

parsedF1L <- parseFile1(file1)

parsedF2L <- parseFile2(file2, SILAC)

parsedF1 <- parsedF1L$parsedF1
minCol <- parsedF1L$minCol
maxCol <- parsedF1L$maxCol
minFrac <- parsedF1L$minFrac
maxFrac <- parsedF1L$maxFrac

alignPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["align"]], sep="")
parsedAlign <- parseAlign(alignPath)

xCalPath <- alignPath <- paste(parsedConstants[["xmldir"]], parsedConstants[["xcal"]], sep="")
parsedxCal <- parsexCal(xCalPath)

trainingSetPath <- paste(scriptDir, trainingSet, sep="/")
svmObject <- svmObjectm(trainingSetPath)

print(paste("Number of cores detected: ", cores, sep=""))

# Temporary holder, used to store peaks detected
peak <- matrix(nrow=parsedF2L$massLength, ncol=20)
colnames(peak) <- c("alRT", "fraction", "peakRTL", "peakRT", "peakRTR", "peakArea",
	"peakWidth", "metric", "sn", "kurtosis", "skew", "derivativeCount", "ms2Call",
	"zeroCount", "maxIntensity", "rtWindow", "peakSig", "peakSharp", "TPASR", "FWHM") #jmb

peak[,"alRT"] <- 0
peak[,"fraction"] <- if(SILAC) -1 else 0
peak[,"peakRTL"] <- 0
peak[,"peakRT"] <- 0
peak[,"peakRTR"] <- 0
peak[,"peakArea"] <- ""
peak[,"peakWidth"] <- 0
peak[,"metric"] <- "no score"
peak[,"sn"] <- ""
peak[,"kurtosis"] <- ""
peak[,"skew"] <- ""
peak[,"derivativeCount"] <-""
peak[,"ms2Call"] <- ""
peak[,"zeroCount"] <- ""
peak[,"maxIntensity"] <- ""
peak[,"rtWindow"] <- ""
peak[,"peakSig"] <- ""
peak[,"peakSharp"] <- ""
peak[,"TPASR"] <- ""
peak[,"FWHM"] <- ""

tempFolder <- parsedConstants[["tempdir"]]

if (as.numeric(cores) < 3) {
    print("Running 1 job.")
    parallel <- FALSE

} else {
    jobs <- cores - 1
    print(paste("Running ", jobs, " jobs in parallel.", sep=""))
    parallel <- TRUE

    cl <- makeCluster(mc <- getOption("cl.cores", jobs), outfile="")

    clusterExport(cl=cl, varlist=c("xcmsScriptPath", "config", "rep", "parsedConstants", "SILAC", 
        "parsedF1", "minCol", "maxCol", "minFrac", "maxFrac", "parsedF2L", "parsedAlign", 
        "parsedxCal", "svmObject", "findEmptyPeaks", "alignRT", "peak", "jobs", "tempFolder", 
        "getRepNumberOfThisRT"))

    invisible(clusterEvalQ(cl, {
        source(xcmsScriptPath)
        suppressWarnings(suppressPackageStartupMessages({
            library(tcltk)
            library(xcms)
            library(e1071)
            library(nortest)
            }))
        }))
}

res <- apply(parsedF1, 1, function(row) {
    column <<- as.numeric(row["column"])
    fraction <<- as.numeric(row["fraction"])
    path <<- as.character(row["path"])
    pathConv <- gsub("\\\\", "/", path)
    fileName <- basename(pathConv)
    fileNameNoExt <- sub("^([^.]*).*", "\\1", fileName)
    mzXMLBase <- paste(fileNameNoExt, ".mzXML", sep="")
    mzXMLPath <<- paste(dirname(path), mzXMLBase, sep="/")
    print(mzXMLPath)
    if (file.exists(mzXMLPath)) {
        file.copy(mzXMLPath, tempFolder)
		print(paste(tempFolder, mzXMLPath, sep=""))
        if (file.exists(paste(tempFolder, mzXMLPath, sep=""))) {
            print("file is moved to right place")
        }
    } else if (file.exists(path)) {
        file.copy(path, tempFolder)
        rawfilePath <- paste(tempFolder, fileName, sep="/")
        #jmb: msconvert call below will be deprecated when mzxml files are available, or will only be run when they're not found
        msconv <- NULL
        if (file.exists(parsedConstants[["msconv"]])) {
            msconv <- parsedConstants[["msconv"]]
        } else if (file.exists(parsedConstants[["msconv64"]])) {
            msconv <- parsedConstants[["msconv64"]]
        }
        if (!is.null(msconv)) {
            cmd <- paste(shQuote(msconv), shQuote(rawfilePath), "-o", shQuote(tempFolder), 
                "--mzXML", "--filter", '"peakPicking true 1-"', sep=" ")
            print(cmd)
            system(cmd, show.output.on.console = TRUE)
        } else {
            stop("msconvert is missing on host computer. Please check ProteoWizard installation.")
        }
    }

    if (file.exists(paste(tempFolder, mzXMLBase, sep="/"))) {
        print(paste("Working on column", column))
        print(paste(tempFolder, mzXMLBase, sep="/"))
        xcmsRawO <<- convert(paste(tempFolder, mzXMLBase, sep="/"))
		print(xcmsRawO)

        if (parallel) {
            clusterExport(cl=cl, varlist=c("xcmsRawO", "column"))

            remove(xcmsRawO)
            gc()

            clusterEvalQ(cl, {
                pb <- tkProgressBar(title=paste("Process:", Sys.getpid()), label=paste("Column:", column), 
                    min=0, max=ceiling(parsedF2L$npeptides/jobs), initial=0, width=300)
                cpeptides <- 0
                })

            # resCol <- parRapply(cl, parsedF2L$parsedF2, function(row) {
            #     cpeptides <<- cpeptides + 1
            #     setTkProgressBar(pb, cpeptides)
            #     return(findEmptyPeaks(row, column, xcmsRawO, peak, parsedF1, 
            #     minCol, maxCol, minFrac, maxFrac, 
            #     as.numeric(parsedConstants[["ppm"]]), as.numeric(rep), 
            #     parsedAlign, parsedxCal, svmObject, SILAC))
            #     })

            # Could not get load balancing working wiht parLapplyLB, had to use clusterApplyLB
            resCol <- clusterApplyLB(cl, 1:parsedF2L$npeptides, function(i) {
                row <<- parsedF2L$parsedF2[i,]
                cpeptides <<- cpeptides + 1
                setTkProgressBar(pb, cpeptides)
                return(findEmptyPeaks(row, column, xcmsRawO, peak, parsedF1, 
                minCol, maxCol, minFrac, maxFrac, 
                as.numeric(parsedConstants[["ppm"]]), as.numeric(rep), 
                parsedAlign, parsedxCal, svmObject, SILAC))
                })

            clusterEvalQ(cl, {
                close(pb)
                remove(xcmsRawO)
                remove(column)
                gc()
                })

        } else {
            pb <- tkProgressBar(title=paste("Process:", Sys.getpid()), label=paste("Column:", column), 
                min=0, max=ceiling(parsedF2L$npeptides), initial=0, width=300)
            cpeptides <- 0

            resCol <- apply(parsedF2L$parsedF2, 1, function(row) {
                cpeptides <<- cpeptides + 1
                setTkProgressBar(pb, cpeptides)
                return(findEmptyPeaks(row, column, xcmsRawO, peak, parsedF1, 
                minCol, maxCol, minFrac, maxFrac, 
                as.numeric(parsedConstants[["ppm"]]), as.numeric(rep), 
                parsedAlign, parsedxCal, svmObject, SILAC))
                })

            close(pb)
            remove(xcmsRawO)
            gc()
        }

        return(do.call(rbind, resCol))
    } else {
        print(paste("Missing raw file:", path))
        return(matrix(list(peak), nrow=parsedF2L$npeptides))
    }
    })

if (parallel) {
    stopCluster(cl) 
}

resF <- lapply(1:parsedF2L$columnLength, function(c) matrix(list(peak), nrow=parsedF2L$npeptides))
xcmsRawCol <- unlist(parsedF1[,"column"])
for (i in 1:length(xcmsRawCol)) {
    resF[xcmsRawCol[i]] <- res[i]
}
resFM <- do.call(cbind, resF)

writeOutput(parsedF2L$parsedF2, resFM, parsedF2L$massLength, parsedF2L$columnLength, parsedF2L$npeptides, outputPath, SILAC)

proc.time() - ptm