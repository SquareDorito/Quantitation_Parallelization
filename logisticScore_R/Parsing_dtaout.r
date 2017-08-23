### R code for parsing in dta and out files 
### To produce computation matrix

curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/Anthony/final_combinations/Rparsing/"

## sets R's current working directory
setwd(curdir)

outfiles <- list.files(pattern=".out")
dtafiles <- list.files(pattern=".dta")

# number of dta files to scan
totalfiles <- length(outfiles)
print(noquote(paste("Processing ",totalfiles, " files")))

# variables to find 
p1 <- vector(mode="numeric", length=totalfiles)
scan1 <- vector(mode="numeric", length=totalfiles)
cs <- vector(mode="numeric", length=totalfiles)
emass <- vector(mode="numeric", length=totalfiles)
 

# generate variables for logistic regression
for (i in 1:totalfiles) {
    currentFile= "BSA_HRP2_LTMSLTMSMS_121405_154555.8055.8055.1.out"
    if (outfiles[i] == currentFile) {                               
        currentFile=outfiles[i]
        print(currentFile)
        #print(outfiles[i])
        ##parses out scan number and charge state            
        p1 <- strsplit(outfiles[i], split=".", fixed=T)         
        scan1[i] <- p1[1]
        cs[i] <- p1[4]
        
        #print(scan1)
        #print(cs)
        
        ## scanning in *.out files
        #out1 <- scan(outfiles[i], what="complex", strip.white=TRUE)
        out1 <- file(description=outfiles[i], open="rt")
        out2 <- readLines(out1)
        
        outtable <- matrix(nrow=10, ncol=8)
        hit1<-matrix(-1,nrow=totalfiles,ncol=7)
        hit2<-matrix(-1,nrow=totalfiles,ncol=7)
        ## reading through each line of every *.out file
        for (j in 1:length(out2)) {       
             str1 <- out2[j]
             print (str1)
             str2 <- strsplit(str1, split=' ',fixed=T)
             #str2<-c(str2)
             #print(str1)
             #print(length(str1))
             #print (j)
             #print(length(str2)) 
             #print(str2)
             

                    
        } 
        
        # closing file connection
        close(out1)
#        temp<-read.delim(currentFile)
#        hit1<-temp[13,1]   #data for hit number 1
#        hit2<-temp[14,1]   #data for hit number 2
#        hit1Stat=array(dim=9)
        
                  
         
    }
}
