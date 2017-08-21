options(warn=-1)
#!/usr/bin/env Rscript
# @aly

########################################################################################
# To make this script run, one must provide three arguments in the following order:
#  1. Raw file path (mzXML must exist in same directory with same basename)
#  2. DTA files folder
#  3. Out files folder
#  4. Peak (output) folder
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
ppm <<- 10

########################################################################################
# Input: OUT file path
# Output: parsed information from OUT file
# 
# This method will parse an outfile and grab information that will be used later on.
#  Specifically it will look for the scan number (from file name). It will also look for 
#  the mass from the line which contains a "1.   ". These are then returned as a list. 
########################################################################################
GetInfoFromOutfile <- function(outfilepath) {
    fNameSplit <- tail(unlist(strsplit(outfilepath, "]")), n=2)
    scan <- as.integer(fNameSplit[1])
    outmass <- NA
    outFile <- file(outfilepath, open="r")
    lines <- readLines(outFile)
    close(outFile)
    for (line in lines) {
        if (grepl("1.   ", line)) {
            s <- gsub("^\\s+|\\s+$", "", line)
            s2 <- unlist(strsplit(s, " +"))
            outmass <- as.numeric(s2[6])
            break
        }
    }
    list(scan=scan, outmass=outmass)
}

########################################################################################
# Input: DTA file path
# Output: parsed information from DTA file
# 
# This method will parse a DTA file and extract the mass and charge from it (first 
#  line).
########################################################################################
GetInfoFromDTAfile <- function(dtafilepath) {
    dtaFile <- file(dtafilepath, open="r")
    lines <- readLines(dtaFile)
    close(dtaFile)
    s <- unlist(strsplit(lines[1], " +"))
    mass <- as.numeric(s[1])
    charge <- as.integer(s[2])
    list(mass=mass, charge=charge)
}

########################################################################################
# Input: xcmsRaw object, output file name, DTA folder path, OUT folder path, peak 
#  (output) folder path
# Output: Nothing
# 
# This method will run the information found using the above methods, and run them 
#  through peakcalc. If a peak is detected, it will write the information to a text file
#  under the name of the out file and dta file (both those files should have the same
#  name). If a peak is not detected, a file is also generated, but with "ND" in all 
#  spots beside the SN, which is set to 1.
# Output file format: Peak intensity, SN, peak width (seconds?), peak apex rt (minutes),
#  all on a new line.
########################################################################################
auto <- function(xcmsRaw, file, dtaFolder, outFolder, peakFolder) {
    outPath <- paste(outFolder, file, sep="/")
    out <- GetInfoFromOutfile(outPath)
    dtaPath <- paste(dtaFolder, file, sep="/")
    dta <- GetInfoFromDTAfile(dtaPath)

    thismass <- out$outmass

    if (thismass == 0) {
        thismass <- (dta$mass - 1.007272) / dta$charge + 1.007272
    } else {
        thismass <- (thismass - 1.007272) / dta$charge + 1.007272
    }

    rt <- acqNumToRt(xcmsRaw, as.numeric(out$scan))

    p <- peakCalc(xcmsRaw, thismass, rt, ppm, NULL, FALSE, TRUE)

    fileOut <- file(paste(peakFolder, basename(file), sep="/"))

    if (is.na(p[["peaks"]][[1]]["peakArea"]) || p[["peaks"]][[1]]["peakArea"] == 0) {
        writeLines(paste("ND", "1", "ND", "ND", sep="\n"), fileOut)
    } else {
        writeLines(paste(p[["peaks"]][[1]]["peakArea"], p[["peaks"]][[1]]["sn"], p[["peaks"]][[1]]["peakWidth"], as.numeric(p[["peaks"]][[1]]["peakApexRT"])/60, sep="\n"), fileOut)
    }
    
    close(fileOut)
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
#   1. Raw file path
#   2. DTA file path (input)
#   3. OUT file path (input)
#   4. peak file path (output)
#   5. msconvert.exe from ProteoWizard path
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

args <- commandArgs(trailingOnly = TRUE)

rawFile <- args[1]
dtaFolder <- args[2]
outFolder <- args[3]
peakFolder <- args[4]
msconv <- args[5]

# dtaFolder <- "E:/Temp/dta"
# outFolder <- "E:/Temp/out"
# peakFolder <- "E:/Temp/peak"

dir.create(file.path(peakFolder), showWarnings = FALSE)
removeOldFiles(peakFolder)

allFiles <- list.files(dtaFolder, include.dirs = FALSE)

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

pb <- tkProgressBar(min=1, max=lengthRows, style=3)

if (as.numeric(cores) < 3) {
    print("Running 1 job.")
    invisible(lapply(allFiles, function(names) {
        write(0, file=f, append=TRUE)
        num <- length(readLines(f, warn = FALSE))
        if (num %% interval == 0) {
            setTkProgressBar(pb, num)
        }
        auto(xcmsRawO, names, dtaFolder, outFolder, peakFolder)
        }))
    # for (names in allFiles) {
    #     write(0, file=f, append=TRUE)
    #     num <- length(readLines(f, warn = FALSE))
    #     if (num %% interval == 0) {
    #         setTxtProgressBar(pb, num)
    #     }
    #     auto(xcmsRawO, names, dtaFolder, outFolder, peakFolder)
    # }
} else {
    jobs <- cores-1
    print(paste("Running ", jobs, " jobs in parallel.", sep=""))
    cl <- makeCluster(mc <- getOption("cl.cores", jobs), outfile="")

    progressdir <- paste(dirname(rawFile), progressFolder, sep="/")
    dir.create(progressdir, showWarnings = FALSE)
    removeOldFiles(progressdir)

    clusterExport(cl=cl, varlist=c("xcmsScriptPath", "xcmsRawO", "allFiles", "dtaFolder", 
        "outFolder", "peakFolder", "f", "lengthRows", "interval", "ppm", 
        "auto", "GetInfoFromOutfile", "GetInfoFromDTAfile", "pb", "progressdir"))

    invisible(clusterEvalQ(cl, {
        source(xcmsScriptPath)
        suppressPackageStartupMessages({
            library(xcms)
            })
        }))

    invisible(clusterApplyLB(cl, allFiles, function(name) {
      fNameSplit <- tail(unlist(strsplit(name, "]")), n=2)
      scan <- as.integer(fNameSplit[1])
      file.create(paste(progressdir, scan, sep=""))
      nFiles <- length(list.files(progressdir, include.dirs = FALSE))
      if (nFiles %% interval == 0) {
        setTkProgressBar(pb, nFiles)
      }
      auto(xcmsRaw=xcmsRawO, file=name, dtaFolder=dtaFolder, outFolder=outFolder, peakFolder=peakFolder)
    }))

    stopCluster(cl)
    removeOldFiles(progressdir)
    gc()
}

proc.time() - ptm