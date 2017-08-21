#!/usr/bin/env Rscript
# @aly

########################################################################################
# To make this script run, one must provide three arguments in the following order:
#  1. CFG file path
#  2. Raw file path (mzXML must exist in same directory with same basename)
#  3. Out files folder
#  4. silac (output) folder
#  5. msconvert.exe from ProteoWizard path
# Also the following packages must be installed:
#  parallel (current versions of R should have it bundled)
#  xcms (from bioconductor)
########################################################################################

suppressPackageStartupMessages({
    library(parallel)
    library(xcms)
    })

########################################################################################
# Names of xcms methods script, msconvert (path), and ppm. Location of them should be in 
#  same folder as script.
########################################################################################
xcmsScript <<- "xcms2.r"
progressFolder <<- "progress/"
# msconv <<- "C:/Program Files (x86)/ProteoWizard/msconvert.exe"
ppm <<- 10 # This does not change the mass window in all cases, only if one is not 
           #  provided in cfg file this is the value that is defaulted to.
aafile <<- "AA_Carbons.txt"

########################################################################################
# Input: CFG file path
# Output: parsed constants
# 
# This method will parse the CFG file as a text file, line by line, and split each line
#  on the "=" symbol. Depending on the text before it, a variable of a similar name will
#  be set. A list of these variables will be returned. For mass window (PPM), if a mass 
#  window is not provided in the file, it will default to 10.
########################################################################################
LoadConfigFromFile <- function(f) {
    label <- list()
    mod <- NA
    output <- NA
    mass_window <- NA
    cfgFile <- file(f, open="r")
    lines <- readLines(cfgFile)
    close(cfgFile)
    labelCounter <- 0
    for (line in lines) {
        split <- unlist(strsplit(line, "="))
        if (split[1] == "LABEL") {
            labelCounter <- labelCounter + 1
            label[[labelCounter]] <- split[2]
        }
        else if (split[1] == "MOD") {
            mod <- as.integer(split[2])
        }
        else if (split[1] == "OUTPUT") {
            output <- split[2]
        }
        else if (split[1] == "MASS_WINDOW") {
            mass_window <- as.integer(split[2])
        }
    }
    if (length(label) == 0) {
        print("No labels in config file.")
    }
    if (is.na(mod)) {
        print("No mod number in config file.")
    }
    if (is.na(output)) {
        print("No output type in config file.")
    }
    if (is.na(mass_window)) {
        mass_window <- ppm
    }
    list(label=label, mod=mod, output=output, mass_window=mass_window)
}

########################################################################################
# Input: OUT file path
# Output: parsed information from OUT file
# 
# This method will parse an OUT file and grab information that will be used later on.
#  Specifically it will look for the scan number (from file name), charge(from file 
#  name), modstring, mass, and peptide sequence. For modstring, it looks for a line that
#  contains "Enzyme: " and grabs the three sets of strings before it, the first two 
#  having parenthesis around it, while the third one is a simple "C=SOME_DOUBLE".
#  Peptide sequence is the amino acid sequence between the periods on the line which 
#  contains a "1.   ", and the mass is also extracted form this line. These are then 
#  returned as a list. 
########################################################################################
GetInfoFromOutfile <- function(outfilepath) {
    fNameSplit <- tail(unlist(strsplit(outfilepath, "]")), n=2)
    scan <- as.integer(fNameSplit[1])
    charge <- as.integer(unlist(strsplit(fNameSplit[2], "\\."))[1])
    modstring <- NA
    outmass <- NA
    pepseq <- NA
    outFile <- file(outfilepath, open="r")
    lines <- readLines(outFile)
    close(outFile)
    for (line in lines) {
        if (grepl("Enzyme:", line)) {
            index <- regexpr("Enzyme:", line, fixed=T)[1]
            ms <- substr(line, 1, index-1)
            ms2 <- gsub("^\\s+|\\s+$", "", ms)
            ms3 <- gsub("\\)", "\\) ", ms2)
            modstring <- unlist(strsplit(ms3, "  +"))
        }
        if (grepl("1.   ", line)) {
            s <- gsub("^\\s+|\\s+$", "", line)
            s2 <- unlist(strsplit(s, " +"))
            outmass <- as.numeric(s2[6])
            pepseq <- unlist(strsplit(tail(s2, n=1), "\\."))[2]
            break
        }
    }
    list(scan=scan, charge=charge, modstring=modstring, outmass=outmass, pepseq=pepseq)
}

########################################################################################
# Input: Labels, mod from CFG, modstring, mass, peptide sequence (last three from OUT 
#  file)
# Output: mass list (2 mass calcualted)
# 
# This method will use the labels and look in the peptide sequence for matches. From 
#  there and the number of matches, it will calculate two masses.
# This is based on the visual basic AQSILAC code.
########################################################################################
InitMassList <- function(isotopelist, isotopesearch, modstring, outmass, pepseq) {
    masslist <- list()
    isotopeL <- length(isotopelist)
    letter <- unlist(strsplit(pepseq, ""))
    letterL <- length(letter)
    aa <- lapply(isotopelist, function(f) unlist(lapply(f, function(g) unlist(strsplit(g, ""))[unlist(gregexpr("\\[", g)[1])+1])))
    dm <- lapply(isotopelist, function(f) unlist(lapply(f, function(g) as.numeric(unlist(strsplit(g, "\\]| |\\|"))[seq(from=2, to=3*length(isotopelist), by=3)]))))
    basemass <- 0
    if (isotopesearch == 0) {
        nomatch <- FALSE
        for (i in 1:isotopeL) {
            tmpaa <- unlist(lapply(1:length(aa[[i]]), function(u) ifelse(CheckDiffMods(aa[[i]][u], dm[[i]][u], modstring) == "x", 
                NA, paste(aa[[i]][u], CheckDiffMods(aa[[i]][u], dm[[i]][u], modstring), sep=""))))
            if (any(is.na(tmpaa))) {
                next
            }
            basemass <- 0
            for (j in 1:letterL) {
                w1 <- FALSE
                if (letter[j] %in% aa[[i]]) {
                    n <- 1
                    if (j < letterL) {
                        if (grepl(letter[j+1], "ARNDCEQGHILKMFPSTYWV")) {
                            n <- 1
                        }
                        else {
                            n <- 2
                        }
                    }
                    for (k in 1:length(tmpaa)) {
                        if (tmpaa[k] == substr(pepseq, j, j+n-1)) {
                            basemass <- basemass - dm[[i]][k]
                            w1 <- TRUE
                            break
                        }
                    }
                    if (w1) {
                        nomatch <- FALSE
                        next
                    }
                    else {
                        nomatch <- TRUE
                        break
                    }
                }
            }
            if (nomatch) {
                next
            }
            else {
                basemass <- outmass + basemass
                for (l in 1:isotopeL) {
                    masslist[[l]] <- basemass
                    for (m in 1:length(aa[[l]])) {
                        for (n in 1:letterL) {
                            if (letter[n] == aa[[l]][m]) {
                                masslist[[l]] <- masslist[[l]] + dm[[l]][m]
                            }
                        }
                    }
                }
                return(list(SILAC=TRUE, masslist=masslist))
            }
        }
    }
    else {
        for (i in 1:isotopeL) {
            masslist[[i]] <- 0
            for (j in 1:length(aa[[i]])) {
                for (k in 1:letterL) {
                    if (letter[k] == aa[[i]][j]) {
                        masslist[[i]] <- masslist[[i]] + dm[[i]][j]
                    }
                }
            }
        }
        basemass <- outmass - masslist[[isotopesearch-1]]
        for (l in 1:isotopeL) {
            masslist[[l]] <- masslist[[l]] + basemass
        }
        return(list(SILAC=TRUE, masslist=masslist))
    }
    return(list(SILAC=FALSE))
}

########################################################################################
# Input: amino acid of label, number associated with label (number in label), modstring
# Output: amino acid with addition (match, no match, etc.)
# 
# This method will first check the value in the label. If it meets the threshold, then
#  it will look for the index where the amino acid falls on the modstring, and check
#  the difference. 
# This method is based on the visual basic AQSILAC.
########################################################################################
CheckDiffMods <- function(a, m, modstring) {
    pos <- lapply(modstring, function(mod) regexpr(a, mod)[1])
    if (abs(m) < 0.2) {
        return("")
    }
    indices <- which(pos != -1)
    if (length(indices) > 0) {
        for (i in indices) {
            num <- as.numeric(gsub("\\)", "", tail(unlist(strsplit(modstring[i], " ")), n=1)))
            if (abs(m - num) < 0.2) {
                return(unlist(strsplit(modstring[i], ""))[regexpr(" ", modstring[i])-1])
            }
        }
    }
    return("x")
}

# SetMassShift <- function(masslist, shiftmass) {
#     for (m in masslist) {
#         m <- m + shiftmass
#     }
#     return(masslist)
# }

# DeIsotopeIt <- function(aapath, isotopelist, pepseq) {
#     # NOT TESTED OR IMPLEMENTED YET
#     cfgFile <- file(aapath, open="r")
#     lines <- readLines(cfgFile)
#     close(cfgFile)
#     aa <- list()
#     cc <- list()
#     for (i in 1:length(lines)) {
#         line <- gsub("^\\s+|\\s+$", "", lines[i])
#         aa[[i]] <- substr(line, 1, 1)
#         cc[[i]] <- substr(line, 2, nchar(line))
#     }

#     il <- isotopelist

#     totalCarbon <- 0
#     j <- 0

#     for (amino in unlist(strsplit(pepseq, ""))) {
#         j <- 0
#         while (j <= length(aa)) {
#             if (amino == aa[[j]]) {
#                 totalCarbon <- totalCarbon + cc[[j]]
#                 break
#             }
#             j <- j + 1
#         }
#     }

# }

########################################################################################
# Input: xcmsRaw object, OUT file path, output file path (saving), parsed CFG file
# Output: Nothing
# 
# This method will run the information found using the above methods, and run them 
#  through peakcalc. Depending on the parameters, it will save the output to one text 
#  or multiple text files under the same name as the input OUT file, but in a different
#  folder. 
# Single text file output not used
# Multiple file output format: 
# First mass peak intensity/SN/mass/stdev peak apex rt/stdev peak shape; Second mass 
#  peak...
########################################################################################
auto <- function(xcmsRaw, outfilepath, outputpath, parsedCfg) {
    out <- GetInfoFromOutfile(outfilepath)
    init <- InitMassList(parsedCfg$label, parsedCfg$mod, out$modstring, out$outmass, out$pepseq)
    masslist <- (unlist(init$masslist)-1.007272) / out$charge + 1.007272
    if (init$SILAC) {
        # if (USE_2ND_ISOTOPIC_PEAK) {
        #     masslist <- SetMassShift(masslist, 1.003355)
        # }
        # if (DEISOTOPE) {

        # }
        if (parsedCfg$output == "txt") {
            rt <- acqNumToRt(xcmsRaw, out$scan)
            p <- peakCalc(xcmsRaw, masslist, rt, parsedCfg$mass_window, NULL, FALSE, FALSE)

            if (substr(basename(outputpath), nchar(basename(outputpath)) - 3, nchar(basename(outputpath))) == ".txt") {
                tmp <- lapply(1:length(masslist), function(y) ifelse(is.na(p[["peaks"]][[y]][["peakArea"]]), 
                    paste("", "",  sep="/"), paste(p[["peaks"]][[y]][["peakArea"]], p[["peaks"]][[y]][["sn"]], sep="/")))
                output <- paste(tmp, collapse = '; ')
                write(output, file=outputpath, append=TRUE)
            }
            else {
                tmp <- lapply(1:length(masslist), function(y) ifelse(is.na(p[["peaks"]][[y]][["peakArea"]]), 
                    paste("", "",  masslist[y], p[["stDevPeakApex"]], p[["stDevPeakShape"]], sep="/"), 
                    paste(p[["peaks"]][[y]][["peakArea"]], p[["peaks"]][[y]][["sn"]], masslist[y], 
                        p[["stDevPeakApex"]], p[["stDevPeakShape"]], sep="/")))
                output <- paste(tmp, collapse = '; ')
                fileOut <- file(paste(outputpath, basename(outfilepath), sep="/"))
                writeLines(output, fileOut)
                close(fileOut)
            }
        }
    }
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
#  Set AA_Carbon.txt path * Not implemented/working
#  Parse arguments:
#   1. CFG file path
#   2. Raw file path
#   3. OUT file path (input)
#   4. silac file path (output)
#   5. msconvert.exe from ProteoWizard path
#  Parse CFG file
#  Check output folder and remove files if there are files in there
#  Get list of input files.
#  Look for mzXML file and parse it into xcmsRaw object
#  Detect number of cores, and run auto on every file in input folder
# 
# Progress bar design: There is a temporary file that for each run of auto
#  (for each row in file 2/peptide), it will write a "0", and count the number of them.
#  Based on the number of them, the progress bar will be updated.
# 
# Multiprocessing design: This uses the built in parallel package bundled with R 2.14+.
#  The current design uses parts of SNOW, and creates a cluster. The parsed objects are 
#  exported to each of the processes in the cluster. Then auto will run in the 
#  same manner as with one process, but with a parallel version of lapply. Also xcmsRaw
#  object from the master process is removed before quantitation and the garbage 
#  collector is invoked to remove it, freeing up memory early on once the xcmsRaw object
#  is copied over to the processes in the cluster. A progress bar is implemented in that
#  it will create blank files based on the scan number and count them.
########################################################################################
ptm <- proc.time()

scriptDir <- currentDir()
xcmsScriptPath <- paste(sep="/", scriptDir, xcmsScript)
cat(paste("Sourcing", xcmsScriptPath, "from", xcmsScriptPath, "\n", sep=" "))
source(xcmsScriptPath)

aapath <- paste(sep="/", scriptDir, aafile)

args <- commandArgs(trailingOnly = TRUE)

configFile <- args[1]
rawFile <- args[2]
inputFolder <- args[3]
outputFolder <- args[4]
msconv <- args[5]

parsedCfg <- LoadConfigFromFile(configFile)

dir.create(file.path(outputFolder), showWarnings = FALSE)
removeOldFiles(outputFolder)

allFiles <- list.files(inputFolder, include.dirs = FALSE)

rawFileDir <- dirname(rawFile)
# cmd <- paste(shQuote(msconv), shQuote(rawFile), "-o", shQuote(rawFileDir), 
#     "--mzXML", "--filter", '"peakPicking true 1-"', sep=" ")
# print(cmd)
# system(cmd, show.output.on.console = TRUE)
fileNameNoExt <- sub("^([^.]*).*", "\\1", basename(rawFile))
# xcmsRawO <- convert(paste(rawFileDir, paste(fileNameNoExt, ".mzXML", sep=""), sep="/"))
xcmsRawO <- convert(paste(rawFileDir, paste(fileNameNoExt, ".mzXML", sep=""), sep="/"))
# file.remove(paste(rawFileDir, paste(fileNameNoExt, ".mzXML", sep=""), sep="/"))

cores <- detectCores()
#cores <- 1
print(paste("Number of cores detected: ", cores, sep=""))

f <- tempfile()
lengthRows <- length(allFiles)
interval <- ceiling(lengthRows/100)

pb <- txtProgressBar(min=1, max=lengthRows, style=3)

if (as.numeric(cores) < 3) {
    print("Running 1 job.")
    invisible(lapply(allFiles, function(names) {
        write(0, file=f, append=TRUE)
        num <- length(readLines(f, warn = FALSE))
        if (num %% interval == 0) {
            setTxtProgressBar(pb, num)
        }
        auto(xcmsRawO, paste(inputFolder, names, sep="/"), outputFolder, parsedCfg)
        }))
    # for (names in allFiles) {
    #     write(0, file=f, append=TRUE)
    #     num <- length(readLines(f, warn = FALSE))
    #     if (num %% interval == 0) {
    #         setTxtProgressBar(pb, num)
    #     }
    #     auto(xcmsRawO, paste(inputFolder, names, sep="/"), outputFolder, parsedCfg)
    # }
} else {
    jobs <- cores-1
    print(paste("Running ", jobs, " jobs in parallel.", sep=""))
    cl <- makeCluster(mc <- getOption("cl.cores", jobs), outfile="")

    progressdir <- paste(dirname(rawFile), progressFolder, sep="/")
    dir.create(progressdir, showWarnings = FALSE)
    removeOldFiles(progressdir)

    clusterExport(cl=cl, varlist=c("xcmsScriptPath", "xcmsRawO", "allFiles", "parsedCfg", 
        "inputFolder", "outputFolder", "f", "lengthRows", "interval", 
        "auto", "LoadConfigFromFile", "GetInfoFromOutfile", "InitMassList", "CheckDiffMods", "progressdir"))

    invisible(clusterEvalQ(cl, {
        source(xcmsScriptPath)
        suppressPackageStartupMessages({
            library(xcms)
            })
        }))

    remove(xcmsRawO)
    gc()

    invisible(clusterApplyLB(cl, allFiles, function(name) {
        fNameSplit <- tail(unlist(strsplit(name, "]")), n=2)
        scan <- as.integer(fNameSplit[1])
        file.create(paste(progressdir, scan, sep=""))
        nFiles <- length(list.files(progressdir, include.dirs = FALSE))
        if (nFiles %% interval == 0) {
            #setTxtProgressBar(pb, nFiles)
        }
        auto(xcmsRaw=xcmsRawO, outfilepath=paste(inputFolder, name, sep="/"), outputpath=outputFolder, parsedCfg=parsedCfg)
        }))

    stopCluster(cl)
    removeOldFiles(progressdir)
}

proc.time() - ptm