currentDir <- function() {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- gsub("\\\\", "/", dirname(script.name))
  return(script.basename)
}

parse_dt_out2<-function(curdir){
#######################################
 ### R code for parsing in dta and out files                       
 ### To produce computation matrix            
 
 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/Anthony/final_combinations/Rparsing/trial/"
 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/matlabandR/Uploaded 12-14/Files/AlphaCasein/acasein_FTMSLTMSMS_HRG2_081105_155721"
 #curdir = "//proteome/User Files/old users/Lisa/mcp5 data/ctrl_IgE_mcp5_5e8"

## sets R's current working directory
  #setwd(curdir)
  print(getwd())

  parsefiles <<- list.files(pattern=".txt")       
  dtafiles <<- list.files(pattern=".dta")
  
  # number of dta files to scan
  totalfiles <<- length(parsefiles)
  print(noquote(paste("Processing ",totalfiles, " files")))
  
  # variables to find
  hit1<<-matrix(-1,9,1)#matrix(-1,nrow=totalfiles,ncol=9)
  hit2<<-matrix(-1,6,1)#matrix(-1,nrow=totalfiles,ncol=6)
  hit<<-1                     
  xtemp<<-0
  ytemp<<-0
  xtempv<<-0
  ytempv<<-0
  xyCoord<<-NULL
  vars<<-NULL
  Fname<<-NULL
  varstemp<<-NULL
  
  suppressPackageStartupMessages({
    suppressWarnings(suppressMessages(library(parallel)))
    suppressWarnings(suppressMessages(library(xcms)))
  })
  
  cores<-detectCores()
  print(paste("Number of cores detected: ", cores,sep=""))
  
  if(as.numeric(cores)<3){
    print("Running 1 job. Not enough cores to parallelize.")
    for (i in 1:totalfiles) {
        
        currentFile= parsefiles[i]#"BSA_HRP2_LTMSLTMSMS_121405_154555.8053.8053.1.out" 
        
        print(noquote(paste(currentFile,": ",i," of ",totalfiles," files")))
        
        par1 <- file(description=parsefiles[i], open="rt")#, encoding="utf-16")
        par2 <- readLines(par1,warn=F)
        dta1 <- file(description=dtafiles[i], open="rt")
        dta2 <- readLines(dta1)
        # 		print(par2)
        # 		print(dta2)  
        if(dta2[2]==""){
          print("fixing .dta file")
          dtaind=1:(length(dta2)/2)
          dtaind=2*dtaind
          dtaind=dtaind-1
          dta2=dta2[dtaind]	
        }
        xtempv=NULL
        ytempv=NULL
        for(a in 1:length(dta2)){
          temp=strsplit(dta2[a],split=' ',fixed=T)
          ytemp=temp[[1]][2] #may be an incorrect label
          xtemp=temp[[1]][1] #may be an incorrect label
          xtempv[a]=as.numeric(xtemp)
          ytempv[a]=as.numeric(ytemp)        
        }
        xyCoord=cbind(xtempv,ytempv)
        dMods=par2[4]
        sMods=par2[5]
        charge=as.numeric(par2[2])
        peptemp=strsplit(par2[6],split="\t",fixed=T)
        #print(par2)
        #print(par2[6])
        peptemp=peptemp[[1]][length(peptemp[[1]])-1]
        pep=strsplit(peptemp,split=".",fixed=T)[[1]][2]
        lineHit2=7
        pep2temp=strsplit(par2[lineHit2],split="\t",fixed=T)
        pep2temp=pep2temp[[1]][length(pep2temp[[1]])-1]
        if(length(pep2temp)>0){
          pep2=strsplit(pep2temp,split=".",fixed=T)[[1]][2] 
          #print(peptemp)
          #print(pep) 
          #print(pep2temp)
          #print(pep2)
          
          while((filterseq(pep)==filterseq(pep2))&&(lineHit2<length(par2))){
            lineHit2=lineHit2+1
            pep2temp=strsplit(par2[lineHit2],split="\t",fixed=T)
            pep2temp=pep2temp[[1]][length(pep2temp[[1]])-1]
            pep2=strsplit(pep2temp,split=".",fixed=T)[[1]][2] 
          }
        }	
        firstHitData=strsplit(par2[6],split="\t",fixed=T)[[1]]
        #print(firstHitData)
        hit1[1]=as.numeric(firstHitData[1])# hit1 (m+h)+
        hit1[2]=as.numeric(firstHitData[2])# hit1 deltCn
        hit1[3]=as.numeric(firstHitData[3])# hit1 XCorr
        hit1[4]=as.numeric(firstHitData[4])# hit1 Sp
        hit1[5]=as.numeric(firstHitData[5])# hit1 Ions Numerator
        hit1[6]=as.numeric(firstHitData[6])# hit1 Ions Denominator
        if(firstHitData[6]==""){hit1[6]=as.numeric(firstHitData[7])} 
        hit1[7]=as.numeric(par2[2])# Charge
        hit1[8]=as.numeric(par2[3])# M+H+mass (exp Mass)
        splitPep=strsplit(pep,split="",fixed=T)[[1]]
        hit1[9]=length(splitPep)-sum(match(splitPep,"[",nomatch=0)+match(splitPep,"]",nomatch=0)+match(splitPep,"*",nomatch=0)+match(splitPep,"#",nomatch=0)+match(splitPep,"@",nomatch=0)+match(splitPep,"^",nomatch=0)+match(splitPep,"$",nomatch=0)+match(splitPep,"~",nomatch=0))# hit1 AA (length of Seq)
        
        secondHitData=strsplit(par2[lineHit2],split="\t",fixed=T)[[1]]
        hit2[1]=as.numeric(secondHitData[1])# hit1 (m+h)+
        hit2[2]=as.numeric(secondHitData[2])# hit1 deltCn
        hit2[3]=as.numeric(secondHitData[3])# hit1 XCorr
        hit2[4]=as.numeric(secondHitData[4])# hit1 Sp
        hit2[5]=as.numeric(secondHitData[5])# hit1 Ions Numerator
        hit2[6]=as.numeric(secondHitData[6])# hit1 Ions Denominator 
        if(is.na(secondHitData[6])||secondHitData[6]==""){hit2[6]=as.numeric(secondHitData[7])}
        
        varstemp[1]=as.numeric(par2[1])#(tempScan[[1]][2])#scan --good
        varstemp[2]=hit1[7]  ##tempScan[[1]][4] #charge --good
        varstemp[3]=xyCoord[1,1] #exp mass  ##Needs to be fixed?
        #print(varstemp[3])
        varstemp[4]=hit1[5] # ions num  --good
        varstemp[5]=hit1[6] # ions den  --good
        varstemp[6]=hit1[4] # sp  ##this is the sp for the first hit, the other code uses the sp from the 10th hit-- dont know why.
        varstemp[7]=hit1[3] # xcorr for hit 1  --good
        varstemp[8]=hit2[3] # xcorr for hit 2  --good
        
        ############################## xcorrPrime calculation start
        
        varstemp[9]=10*log(varstemp[7])/(log(2*varstemp[2]*(hit1[9]-1))) #xcorrprime  ##AA in paper AA-1 in code?? --good
        
        ############################ xcorrPrime calculation end
        varstemp[10]=hit2[2]*10 #delta cn hit2  --good
        if(is.na(varstemp[10])){
          varstemp[10]=1		# if no second hit, set dcn2 to 1 arbitrary
        }
        varstemp[11]=hit1[1] # mhmass(1) --good
        #print(varstemp[11])
        varstemp[12]=10*hit1[5]/hit1[6] #ionsratio    ## 100* in paper 10* in code?? --good
        varstemp[13]=hit1[9] #AA (number of amino acids) --good
        
        ############################ delta mass 1 and delta mass 2 calc start
        
        varstemp[14]=(1e6)*abs((xyCoord[1,1]-hit1[1])/xyCoord[1,1]) #dmass1      --good!!
        #print(varstemp[14])
        #bob=(1e6)*abs((xyCoord[1,1]-hit1[i,1])/hit1[i,1])
        #print(bob)
        varstemp[15]=(1e6)*abs((xyCoord[1,1]-hit2[1])/xyCoord[1,1]) #dmass2      --good!!
        if(is.na(varstemp[15])){
          varstemp[15]=1000		# if no second hit, set dmass2 to 1000ppm arbitrary
        }
        ########################## delta mass 1 and delta mass 2 calc end
        
        #varstemp[15]=(hit1[i,3]-hit2[i,3])/hit1[i,3] # dc2 or delta c2
        
        ######################### start calculate number of K and R in sequence   #--good
        numK=length(splitPep)-sum(as.numeric(is.na(pmatch(splitPep,"K",duplicates.ok=TRUE))))
        numR=length(splitPep)-sum(as.numeric(is.na(pmatch(splitPep,"R",duplicates.ok=TRUE))))
        ## there is a cleaner way to do this using logic in the phosSites calc
        varstemp[16]=numK+numR-1 #KR  --good
        #print(varstemp[16])
        ## note that in the matlab code we subract 1 from this value; not sure if that is needed 
        ## here, but i think it is to account for the theoretical terminal K/R
        
        ######################## end calculate number of K and R in sequence
        
        ##################### ps pt py and phos sites calc start        
        # !!!!! #May 29 2008.  Resume correcting below.    
        dm=getDMods(dMods)
        logical1=dm>79
        logical2=dm<=80
        phosSymbol=names(dm[logical1&logical2]) 
        numPS=0
        numPT=0
        numPY=0
        if(length(phosSymbol)>0){ 
          phosIndex=match(splitPep,phosSymbol,nomatch=0)
          phosIndex=as.logical(c(phosIndex[2:length(phosIndex)],phosIndex[1]))
          phosSites=splitPep[phosIndex]
          if(length(phosSites)>0){
            for(j in 1:length(phosSites)){
              if(phosSites[j]=="S"){numPS=numPS+1}
              if(phosSites[j]=="T"){numPT=numPT+1}
              if(phosSites[j]=="Y"){numPY=numPY+1}
            }
          }
        }
        
        
        varstemp[17]=numPS+numPT+numPY
        varstemp[18]=numPS  #--good
        varstemp[19]=numPT  #--good
        varstemp[20]=numPY  #--good(assumed as previous 2 were good)
        
        #################### ps pt py and phos sites calc end
        
        #################### start number mean number fancymean median calc  ## see noise.m
        
        no=length(xyCoord[,2])-1 #subtract the peak for exp mass
        sorted=sort(xyCoord[2:leng(xyCoord),2],T)#sort peaks note:F=ascending order
        #print(sorted)
        normalized=xyCoord[2:leng(xyCoord),2]/max(xyCoord[2:leng(xyCoord),2]) #--need to eliminate 1st entry
        xbarfull=mean(xyCoord[2:leng(xyCoord),2])
        xbarminus=mean(sorted[(varstemp[13]+1):no])
        #print(xbarfull)
        #print(xbarminus)
        varstemp[21]=(xbarfull-xbarminus)/xbarfull #fancymean  --good                
        varstemp[22]=no/varstemp[13] #number    --good
        #print(varstemp[21])
        #print(varstemp[22])
        varstemp[23]=median(normalized) #median --good   
        varstemp[24]=mean(normalized) #mean     --good   
        
        #################### end number mean number fancymean median calc                                                  
        
        #################### bscore & sumscore calc start
        
        ## note that hit1[i,8] is supposed to be what is MHMass(1)
        ## in the matlab code seqprocess.m and varstemp[2] corresponds to
        ## chargeState in that same code
        isoMass=(varstemp[11]+1)/varstemp[2]     ##-good
        MHH3PO4=isoMass-(98/varstemp[2])##MminusH(1)    ##-good
        MHHPO3=isoMass-(80/varstemp[2])          ##-good
        MH2H3PO4=MHH3PO4-(98/varstemp[2])        ##-good
        MH2HPO3=MHHPO3-(80/varstemp[2])          ##-good
        
        # WORKS! #print((isolate(xyCoord[2:length(xyCoord[,1]),])))
        #print(xyCoord)
        loners=isolate(xyCoord)  #-good
        #print(xyCoord)
        #print(loners[[1]])
        #        print(min(varstemp[13],length(loners[[1]][,1])))
        
        top=findTop(min(varstemp[13],length(loners[[1]][,1])),loners[[1]]) #-good
        #print(top)               
        ##RESUME DEBUGGING HERE
        ##the folowing is from searchcp.m 
        if(varstemp[17]==0){
          varstemp[25]=0#sumscore
          varstemp[26]=0#bscore
        }
        else{
          ub=MHH3PO4+.5
          lb=MHH3PO4-.5
          wloss=18/varstemp[2]
          aloss=17/varstemp[2]
          mloss=98/varstemp[2]
          locs=NULL
          locs[[1]]=max(vectorFind((top[,1]<=ub)&(top[,1]>=lb)),0)
          locs[[2]]=max(vectorFind((top[,1]<=ub-wloss)&(top[,1]>=lb-wloss)),0)
          locs[[3]]=max(vectorFind((top[,1]<=ub-aloss)&(top[,1]>=lb-aloss)),0)
          locs[[4]]=max(vectorFind((top[,1]<=ub-(2*wloss))&(top[,1]>=lb-(2*wloss))),0)
          locs[[5]]=max(vectorFind((top[,1]<=ub-(2*aloss))&(top[,1]>=lb-(2*aloss))),0)
          #print(locs)
          found=array(0,5)
          inten=array(0,5)
          
          for(k in 1:5){
            if((locs[[k]])!=0){
              found[k]=1
              peaks=NULL
              #print(locs[[k]])
              for(a in 1:length(locs[[k]])){
                if(locs[[k]][a]!=0){peaks=cbind(peaks,t(as.matrix(top[locs[[k]][a],])))}   ##changed!!                                
              }
              inten[k]=max(peaks[,2])
            }
          }
          
          varstemp[25] = (2*found[1])+(found[2]+found[3])+(.5*sum(found[4:5]))#bscore   --(not sure)--presumed correct
          varstemp[26]=sum(inten)#sumscore                                              --(not sure)--presumed correct
        }
        #################### bscore & sumscore calc end
        
        #################### noa nda nsa toa tda tsa etc... calc start
        pepc=pepcalc(pep,charge,dMods,sMods)
        #print(pepc)
        MminusH=c(MHH3PO4,MH2H3PO4)
        
        us=uscore(top,pepc,MminusH,varstemp[2],varstemp[17],phosSymbol)
        
        
        #################### noa nda nsa toa tda tsa etc... calc end
        if(firstHitData[8]=="R"){reDb=1}
        else{reDb=0}
        
        vars[[i]]=cbind(as.matrix(us),t(as.matrix(varstemp)),reDb)
        Fname[i]=parsefiles[i]
        close(par1)
        close(dta1)
      }
  }else{
    jobs<-cores-1
    print(paste("Running ",jobs," jobs in parallel.",sep=""))
    #print("hello1")
    cl <- makeCluster(mc <- getOption("cl.cores", jobs), outfile="")
    #print("hello2")
    clusterExport(cl=cl, varlist=c("parsefiles","Fname","varstemp","vars",
                                  "hit1","hit2","xtemp","ytemp","xtempv","ytempv","xyCoord",
                                  "hit","totalfiles","dtafiles","getDMods","leng","isolate",
                                  "vectorFind","findTop","pepcalc","uscore","getSMods",
                                  "getavgmass","getmass","filterseq"))
    #print("hello3")
    invisible(clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        suppressWarnings(suppressMessages(library(xcms)))
      })
    }))
    #print("hello4")
    invisible(clusterApplyLB(cl,parsefiles,function(name){
      #print("hello5")
      counter<-match(c(name),parsefiles)
      print(counter,sep="")
      currentFile=name
      par1<-file(description=currentFile,open="rt")
      par2<-readLines(par1,warn=F)
      #print(par2,sep="\n")
      dta1<-file(description=dtafiles[counter],open="rt")
      dta2<-readLines(dta1)
      #print(dta2,sep="\n")
      
      if(dta2[2]==""){
        print("fixing .dta file")
        dtaind=1:(length(dta2)/2)
        dtaind=2*dtaind
        dtaind=dtaind-1
        dta2=dta2[dtaind]	
      }
      xtempv=NULL
      ytempv=NULL
      for(a in 1:length(dta2)){
        temp=strsplit(dta2[a],split=' ',fixed=T)
        ytemp=temp[[1]][2] #may be an incorrect label
        xtemp=temp[[1]][1] #may be an incorrect label
        xtempv[a]=as.numeric(xtemp)
        ytempv[a]=as.numeric(ytemp)        
      }
      xyCoord=cbind(xtempv,ytempv)
      dMods=par2[4]
      sMods=par2[5]
      charge=as.numeric(par2[2])
      peptemp=strsplit(par2[6],split="\t",fixed=T)
      #print(par2)
      #print(par2[6])
      peptemp=peptemp[[1]][length(peptemp[[1]])-1]
      pep=strsplit(peptemp,split=".",fixed=T)[[1]][2]
      lineHit2=7
      pep2temp=strsplit(par2[lineHit2],split="\t",fixed=T)
      pep2temp=pep2temp[[1]][length(pep2temp[[1]])-1]
      if(length(pep2temp)>0){
        pep2=strsplit(pep2temp,split=".",fixed=T)[[1]][2] 
        #print(peptemp)
        #print(pep) 
        #print(pep2temp)
        #print(pep2)
        
        while((filterseq(pep)==filterseq(pep2))&&(lineHit2<length(par2))){
          lineHit2=lineHit2+1
          pep2temp=strsplit(par2[lineHit2],split="\t",fixed=T)
          pep2temp=pep2temp[[1]][length(pep2temp[[1]])-1]
          pep2=strsplit(pep2temp,split=".",fixed=T)[[1]][2] 
        }
      }
      firstHitData=strsplit(par2[6],split="\t",fixed=T)[[1]]
      #print(firstHitData)
      hit1[1]=as.numeric(firstHitData[1])# hit1 (m+h)+
      hit1[2]=as.numeric(firstHitData[2])# hit1 deltCn
      hit1[3]=as.numeric(firstHitData[3])# hit1 XCorr
      hit1[4]=as.numeric(firstHitData[4])# hit1 Sp
      hit1[5]=as.numeric(firstHitData[5])# hit1 Ions Numerator
      hit1[6]=as.numeric(firstHitData[6])# hit1 Ions Denominator
      if(firstHitData[6]==""){hit1[6]=as.numeric(firstHitData[7])} 
      hit1[7]=as.numeric(par2[2])# Charge
      hit1[8]=as.numeric(par2[3])# M+H+mass (exp Mass)
      splitPep=strsplit(pep,split="",fixed=T)[[1]]
      hit1[9]=length(splitPep)-sum(match(splitPep,"[",nomatch=0)+match(splitPep,"]",nomatch=0)+match(splitPep,"*",nomatch=0)+match(splitPep,"#",nomatch=0)+match(splitPep,"@",nomatch=0)+match(splitPep,"^",nomatch=0)+match(splitPep,"$",nomatch=0)+match(splitPep,"~",nomatch=0))# hit1 AA (length of Seq)
      
      secondHitData=strsplit(par2[lineHit2],split="\t",fixed=T)[[1]]
      hit2[1]=as.numeric(secondHitData[1])# hit1 (m+h)+
      hit2[2]=as.numeric(secondHitData[2])# hit1 deltCn
      hit2[3]=as.numeric(secondHitData[3])# hit1 XCorr
      hit2[4]=as.numeric(secondHitData[4])# hit1 Sp
      hit2[5]=as.numeric(secondHitData[5])# hit1 Ions Numerator
      hit2[6]=as.numeric(secondHitData[6])# hit1 Ions Denominator 
      if(is.na(secondHitData[6])||secondHitData[6]==""){hit2[6]=as.numeric(secondHitData[7])}
      
      varstemp[1]=as.numeric(par2[1])#(tempScan[[1]][2])#scan --good
      varstemp[2]=hit1[7]  ##tempScan[[1]][4] #charge --good
      varstemp[3]=xyCoord[1,1] #exp mass  ##Needs to be fixed?
      #print(varstemp[3])
      varstemp[4]=hit1[5] # ions num  --good
      varstemp[5]=hit1[6] # ions den  --good
      varstemp[6]=hit1[4] # sp  ##this is the sp for the first hit, the other code uses the sp from the 10th hit-- dont know why.
      varstemp[7]=hit1[3] # xcorr for hit 1  --good
      varstemp[8]=hit2[3] # xcorr for hit 2  --good
      
      ############################## xcorrPrime calculation start
      
      varstemp[9]=10*log(varstemp[7])/(log(2*varstemp[2]*(hit1[9]-1))) #xcorrprime  ##AA in paper AA-1 in code?? --good
      
      ############################ xcorrPrime calculation end
      varstemp[10]=hit2[2]*10 #delta cn hit2  --good
      if(is.na(varstemp[10])){
        varstemp[10]=1		# if no second hit, set dcn2 to 1 arbitrary
      }
      varstemp[11]=hit1[1] # mhmass(1) --good
      #print(varstemp[11])
      varstemp[12]=10*hit1[5]/hit1[6] #ionsratio    ## 100* in paper 10* in code?? --good
      varstemp[13]=hit1[9] #AA (number of amino acids) --good
      
      ############################ delta mass 1 and delta mass 2 calc start
      
      varstemp[14]=(1e6)*abs((xyCoord[1,1]-hit1[1])/xyCoord[1,1]) #dmass1      --good!!
      #print(varstemp[14])
      #bob=(1e6)*abs((xyCoord[1,1]-hit1[i,1])/hit1[i,1])
      #print(bob)
      varstemp[15]=(1e6)*abs((xyCoord[1,1]-hit2[1])/xyCoord[1,1]) #dmass2      --good!!
      if(is.na(varstemp[15])){
        varstemp[15]=1000		# if no second hit, set dmass2 to 1000ppm arbitrary
      }
      ########################## delta mass 1 and delta mass 2 calc end
      
      #varstemp[15]=(hit1[i,3]-hit2[i,3])/hit1[i,3] # dc2 or delta c2
      
      ######################### start calculate number of K and R in sequence   #--good
      numK=length(splitPep)-sum(as.numeric(is.na(pmatch(splitPep,"K",duplicates.ok=TRUE))))
      numR=length(splitPep)-sum(as.numeric(is.na(pmatch(splitPep,"R",duplicates.ok=TRUE))))
      ## there is a cleaner way to do this using logic in the phosSites calc
      varstemp[16]=numK+numR-1 #KR  --good
      #print(varstemp[16])
      ## note that in the matlab code we subract 1 from this value; not sure if that is needed 
      ## here, but i think it is to account for the theoretical terminal K/R
      
      ######################## end calculate number of K and R in sequence
      
      ##################### ps pt py and phos sites calc start        
      # !!!!! #May 29 2008.  Resume correcting below.    
      dm=getDMods(dMods)
      #print("after d mods",sep="\n")
      logical1=dm>79
      logical2=dm<=80
      phosSymbol=names(dm[logical1&logical2]) 
      numPS=0
      numPT=0
      numPY=0
      if(length(phosSymbol)>0){ 
        phosIndex=match(splitPep,phosSymbol,nomatch=0)
        phosIndex=as.logical(c(phosIndex[2:length(phosIndex)],phosIndex[1]))
        phosSites=splitPep[phosIndex]
        if(length(phosSites)>0){
          for(j in 1:length(phosSites)){
            if(phosSites[j]=="S"){numPS=numPS+1}
            if(phosSites[j]=="T"){numPT=numPT+1}
            if(phosSites[j]=="Y"){numPY=numPY+1}
          }
        }
      }
      
      
      varstemp[17]=numPS+numPT+numPY
      varstemp[18]=numPS  #--good
      varstemp[19]=numPT  #--good
      varstemp[20]=numPY  #--good(assumed as previous 2 were good)
      
      #################### ps pt py and phos sites calc end
      
      #################### start number mean number fancymean median calc  ## see noise.m
      
      no=length(xyCoord[,2])-1 #subtract the peak for exp mass
      sorted=sort(xyCoord[2:leng(xyCoord),2],T)#sort peaks note:F=ascending order
      #print(sorted,sep="\n")
      normalized=xyCoord[2:leng(xyCoord),2]/max(xyCoord[2:leng(xyCoord),2]) #--need to eliminate 1st entry
      xbarfull=mean(xyCoord[2:leng(xyCoord),2])
      xbarminus=mean(sorted[(varstemp[13]+1):no])
      #print(xbarfull)
      #print(xbarminus)
      varstemp[21]=(xbarfull-xbarminus)/xbarfull #fancymean  --good                
      varstemp[22]=no/varstemp[13] #number    --good
      #print(varstemp[21])
      #print(varstemp[22])
      varstemp[23]=median(normalized) #median --good   
      varstemp[24]=mean(normalized) #mean     --good   
      
      #################### end number mean number fancymean median calc                                                  
      
      #################### bscore & sumscore calc start
      
      ## note that hit1[i,8] is supposed to be what is MHMass(1)
      ## in the matlab code seqprocess.m and varstemp[2] corresponds to
      ## chargeState in that same code
      isoMass=(varstemp[11]+1)/varstemp[2]     ##-good
      MHH3PO4=isoMass-(98/varstemp[2])##MminusH(1)    ##-good
      MHHPO3=isoMass-(80/varstemp[2])          ##-good
      MH2H3PO4=MHH3PO4-(98/varstemp[2])        ##-good
      MH2HPO3=MHHPO3-(80/varstemp[2])          ##-good
      
      # WORKS! #print((isolate(xyCoord[2:length(xyCoord[,1]),])))
      #print(xyCoord)
      loners=isolate(xyCoord)  #-good
      #print(xyCoord)
      #print(loners[[1]])
      #        print(min(varstemp[13],length(loners[[1]][,1])))
      
      top=findTop(min(varstemp[13],length(loners[[1]][,1])),loners[[1]]) #-good
      #print(top)               
      ##RESUME DEBUGGING HERE
      ##the folowing is from searchcp.m 
      if(varstemp[17]==0){
        varstemp[25]=0#sumscore
        varstemp[26]=0#bscore
      }
      else{
        ub=MHH3PO4+.5
        lb=MHH3PO4-.5
        wloss=18/varstemp[2]
        aloss=17/varstemp[2]
        mloss=98/varstemp[2]
        locs=NULL
        locs[[1]]=max(vectorFind((top[,1]<=ub)&(top[,1]>=lb)),0)
        locs[[2]]=max(vectorFind((top[,1]<=ub-wloss)&(top[,1]>=lb-wloss)),0)
        locs[[3]]=max(vectorFind((top[,1]<=ub-aloss)&(top[,1]>=lb-aloss)),0)
        locs[[4]]=max(vectorFind((top[,1]<=ub-(2*wloss))&(top[,1]>=lb-(2*wloss))),0)
        locs[[5]]=max(vectorFind((top[,1]<=ub-(2*aloss))&(top[,1]>=lb-(2*aloss))),0)
        #print(locs)
        found=array(0,5)
        inten=array(0,5)
        
        for(k in 1:5){
          if((locs[[k]])!=0){
            found[k]=1
            peaks=NULL
            #print(locs[[k]])
            for(a in 1:length(locs[[k]])){
              if(locs[[k]][a]!=0){peaks=cbind(peaks,t(as.matrix(top[locs[[k]][a],])))}   ##changed!!                                
            }
            inten[k]=max(peaks[,2])
          }
        }
        
        varstemp[25] = (2*found[1])+(found[2]+found[3])+(.5*sum(found[4:5]))#bscore   --(not sure)--presumed correct
        varstemp[26]=sum(inten)#sumscore                                              --(not sure)--presumed correct
      }
      #################### bscore & sumscore calc end
      
      #################### noa nda nsa toa tda tsa etc... calc start
      pepc=pepcalc(pep,charge,dMods,sMods)
      #print(pepc)
      MminusH=c(MHH3PO4,MH2H3PO4)
      
      us=uscore(top,pepc,MminusH,varstemp[2],varstemp[17],phosSymbol)
      
      
      #################### noa nda nsa toa tda tsa etc... calc end
      if(firstHitData[8]=="R"){reDb=1}
      else{reDb=0}
      
      #write.table(cbind(as.matrix(us),t(as.matrix(varstemp)),reDb),file=paste("testing/",counter,"_vars.txt",sep=""),sep="\t")
      #write.table(name,file=paste("testing/",counter,"_fname.txt",sep=""),col.names=FALSE,row.names=FALSE)
      vars[[counter]]<<-cbind(as.matrix(us),t(as.matrix(varstemp)),reDb)
      Fname[counter]<<-name
      write.table(vars[[counter]],file=paste("testing/vars",counter,".txt",sep=""),sep="\t",col.names=FALSE)
      write(Fname[counter],file="testing/fname.txt",append=TRUE)
      #print(name)
      #print("DONE",sep="\n")
      close(par1)
      close(dta1)
    }))
    #todo
  }
  
  # generate variables for logistic regression
  retVars=NULL
  for(i in 1:totalfiles){
    write.table(vars[[i]],file=paste("testing/aftervars",i,".txt",sep=""),sep="\t",col.names=FALSE)
  }
  # tempretvars=NULL
  # for(i in 1:totalfiles){
  #   tempretvars[i]=read.table(paste("testing/vars",i,".txt",sep=""),sep="\t")
  # }
  # write.table(tempretvars[1],file=paste("testing/pvars1.txt",sep=""),sep="\t")
  # fnames<-read.table("testing/fname.txt",sep=",")
  # for(i in 1:totalfiles){
  #   Fname[i]=fnames[i,1]
  #   
  # }
  #write(Fname,file="processedfname.txt",sep=",")
  retVars[[1]]=vars
  retVars[[2]]=Fname
  
  stopCluster(cl) 
  gc()
  
  return(retVars)
}

filterseq<-function(pepseq){
	var=gsub("[[:punct:]]","",pepseq)
	return(var)
}
            
            
            
            
            
            
            
            
            