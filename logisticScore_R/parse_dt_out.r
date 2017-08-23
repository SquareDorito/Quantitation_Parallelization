parse_dt_out<-function(curdir){
#######################################
 ### R code for parsing in dta and out files                       
 ### To produce computation matrix            
 
 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/Anthony/final_combinations/Rparsing/trial/"
 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/matlabandR/Uploaded 12-14/Files/AlphaCasein/acasein_FTMSLTMSMS_HRG2_081105_155721"
 #curdir = "//proteome/User Files/old users/Lisa/mcp5 data/ctrl_IgE_mcp5_5e8"

## sets R's current working directory
#setwd(curdir)

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
hit1<-matrix(-1,nrow=totalfiles,ncol=9)
hit2<-matrix(-1,nrow=totalfiles,ncol=6)
hit=1                     
xtemp=0
ytemp=0
xtempv=0
ytempv=0
xyCoord=NULL
vars=NULL
varstemp=NULL
# generate variables for logistic regression
for (i in 1:totalfiles) {

    currentFile= outfiles[i]#"BSA_HRP2_LTMSLTMSMS_121405_154555.8053.8053.1.out" 
    
        print(noquote(paste(currentFile,": ",i," of ",totalfiles," files")))

        out1 <- file(description=outfiles[i], open="rt")
        out2 <- readLines(out1)
        dta1 <- file(description=dtafiles[i], open="rt")
        dta2 <- readLines(dta1)
        #print(out1)
        #print(out2)
        if(length(out2)>35){ ##if unix style line endings
          theseLines=1:floor(length(out2)/2)  
          thoseLines=1:(length(dta2)/2)
          out2=out2[2*theseLines-1]
          dta2=dta2[2*thoseLines-1]
        }
        #print(out2)
        
        print(out2)
        ## reading through each line of every *.out file
        #for (j in 1:length(out2)) {
        xtempv=NULL
        ytempv=NULL
        for(a in 1:length(dta2)){
            temp=strsplit(dta2[a],split=' ',fixed=T)
            ytemp=temp[[1]][2] #may be an incorrect label
            xtemp=temp[[1]][1] #may be an incorrect label
            xtempv[a]=as.numeric(xtemp)
            ytempv[a]=as.numeric(ytemp)        
        }
        xyCoord[[i]]=cbind(xtempv,ytempv)
        #print(xyCoord[[i]])
        #print(xyCoord[[i]])
        ########
        j=12  
        str1<-out2[j]
        #print(str1)
        mods=outmodsFile(str1)
        print(str1)
        print(mods)
        ########
             hit=NULL 
             j=16
             str1 <- out2[j]
             admit=1
             if(j==16) {
                str2 <- strsplit(str1, split=' ',fixed=T)
                count<-1
                for (k in 1:length(str2[[1]])){
                if(length(str2[[1]])<=1){
                  admit=0
                  #print('here')
                  }
                  else if(str2[[1]][k]!='') {
                       temp<-strsplit(str2[[1]][k], split='/',fixed=T)
                       #print(temp)
                       for(m in 1:length(temp[[1]])){
                           if(temp[[1]][m]!='') {
                              hit[count]<-temp[[1]][m]
                              count<-count+1
                           }
                       }
                    }
                
                }
                #print(hit)
                if(admit==0){
                vars[[i]]=c(outfiles[i],dtafiles[i])
                #print('here_1')
                #close(out1)
                #close(dta1)
                
                }
                else{
                hit1[i,1]<-as.numeric(hit[5]) # hit1 (m+h)+
                hit1[i,2]<-as.numeric(hit[6]) # hit1 deltCn
                hit1[i,3]<-as.numeric(hit[7]) # hit1 xcorr
                hit1[i,4]<-as.numeric(hit[8]) # hit1 sp
                hit1[i,5]<-as.numeric(hit[9]) # hit1 ionsNum
                hit1[i,6]<-as.numeric(hit[10]) # hit1 ionsDen
                #hit1[i,7]<-as.numeric(hit[12]) #hit1 charge
                str2=strsplit(out2[7],split='+',fixed=T)
                str2=strsplit(str2[[1]][4],split='',fixed=T)
                hit1[i,7]<-as.numeric(str2[[1]][1]) #hit1 charge
                str2=strsplit(out2[7],split='=',fixed=T)
                str2=strsplit(str2[[1]][2],split=' ',fixed=T)
                hit1[i,8]<-as.numeric(str2[[1]][2]) #hit1 (m+h)+mass
                lengthTemp=strsplit(hit[length(hit)],split='.',fixed=T)
                posTemp=lengthTemp[[1]][2]
                lengthTemp=strsplit(lengthTemp[[1]][2],split='',fixed=T)
                posTempSep=strsplit(posTemp,split='',fixed=T)                     
                phosIndex=pmatch(posTempSep[[1]],"*",nomatch=0,duplicates.ok=TRUE)
                #print(lengthTemp)
                #print(phosIndex)
                hit1[i,9]=length(lengthTemp[[1]])-sum(phosIndex) #hit1 #AA   ##subtraction accounts for number of * in string
                hit=NULL
                firstSeq2=strsplit(posTemp,"*",fixed=T)
                firstSeq=NULL
                for(z in 1:length(firstSeq2[[1]])){
                   firstSeq=paste(firstSeq,firstSeq2[[1]][z],sep="")
                }
                }
             }
             if(admit==1){admit=filterFile(posTemp,hit1[i,3],hit1[i,7])}
             if(admit==0){
             vars[[i]]=c(outfiles[i],dtafiles[i])
             #print('here_2')
             #close(out1)
             #close(dta1)
             }
             else{
     #        ################################## here we iterate to find next seq diff hit         INCOMPLETE
             flag=T
             while(flag){    ##NEED TO FIX FOR CASE WHERE ALL HITS ARE SAME SEQUENCE!   --DONE!
               j=j+1
               if(j==26){
                  j=17
                  flag=FALSE
                  #print("here1")
                  break
                  #print("here2")
               }
              temp<-strsplit(out2[j], split=' ',fixed=T)
                         #print(temp)
              seqTemp=strsplit(temp[[1]][length(temp[[1]])],'.',fixed=T)
              seqTemp2=strsplit(seqTemp[[1]][2],"*",fixed=T)
              seqTemp=NULL
              for(y in 1:length(seqTemp2[[1]])){
                 seqTemp=paste(seqTemp,seqTemp2[[1]][y],sep="")
              }
              if(seqTemp!=firstSeq){flag=FALSE}
              #if(j==25){
              #  flag=FALSE
              #  j=17
              #}
             }                                                                         
     #        #########################################################################
             #j=17                #doctored for above
             str1 <- out2[j]
             if(!flag) {          # "        "   "
                str2 <- strsplit(str1, split=' ',fixed=T)
                count<-1
                for (k in 1:length(str2[[1]])){
                    if(str2[[1]][k]!='') {
                       temp<-strsplit(str2[[1]][k], split='/',fixed=T)
                       #print(temp)
                       for(m in 1:length(temp[[1]])){
                           if(temp[[1]][m]!='') {
                              hit[count]<-temp[[1]][m]
                              count<-count+1
                           }
                       }
                    }

                }
                
                hit2[i,1]<-as.numeric(hit[5]) # hit2 mhmass
                hit2[i,2]<-as.numeric(hit[6]) # hit2 deltCn
                hit2[i,3]<-as.numeric(hit[7]) # hit2 xcorr
                hit2[i,4]<-as.numeric(hit[8]) # hit2 sp
                hit2[i,5]<-as.numeric(hit[9]) # hit2 ionsNum
                hit2[i,6]<-as.numeric(hit[10]) # hit2 ionsDen
                hit=NULL
             }
             
         #################### here we are going to patch in a new sp, ions num and ions den to 
         #################### be the same as MATLAB-- this may need to be removed to 
         #################### improve results

 #         j=25
 #         count=1
 #         str1=out2[j]
 #         str2 <- strsplit(str1, split=' ',fixed=T)
 #         for (k in 1:length(str2[[1]])){                            
 #         if(length(str2[[1]])<=1){                                  
 #           admit=0                                                  
 #           #print('here')                                           
 #           }                                                        
 #           else if(str2[[1]][k]!='') {                              
 #                temp<-strsplit(str2[[1]][k], split='/',fixed=T)     
 #                #print(temp)                                        
 #                for(m in 1:length(temp[[1]])){                      
 #                    if(temp[[1]][m]!='') {                          
 #                       hit[count]<-temp[[1]][m]                     
 #                       count<-count+1                               
 #                    }                                               
 #                }                                                   
 #             }                                                      
 #                                                                    
 #         }                                                          
 #          hit1[i,4]<-as.numeric(hit[8]) # hit1 sp        
 #          hit1[i,5]<-as.numeric(hit[9]) # hit1 ionsNum   
 #          hit1[i,6]<-as.numeric(hit[10]) # hit1 ionsDen  
          
          #lastSeqLine=strsplit(out2[j],split='/',fixed=T)
          #numLine=strsplit(lastSeqLine[[1]][2],split=' ',fixed=T)
          #denLine=strsplit(lastSeqLine[[1]][3],split=' ',fixed=T)
          #hit1[i,5]=as.numeric(numLine[[1]][length(numLine[[1]])])
          #hit1[i,6]=as.numeric(denLine[[1]][1])
          #spInd=as.numeric(is.na(as.numeric(numLine[[1]])))
          #spInd=vectorFind(as.numeric(spInd==0))
          #hit1[i,4]=as.numeric(numLine[[1]][spInd[length(spInd)-1]]) 
          
         #################### end patch    
          
        #}
          
        # closing file connection
        #print('here_3')
        #close(out1)
        #close(dta1)

        ## hit1 stores data for first hit by row
        ## hit2 stores data for second hit by row
        ## xyCoord stores coordinates of peaks by cell
        ##       i.e. xyCoord[[1]] is the first file and xyCoord[[1]][2,1] is an
        ##            x cood. and xyCoord[[1]][2,2] is a y coord.
        ##            xyCoord[[a]][1,1] is the precursor mass for the a-th file.
        tempScan=strsplit(currentFile,split='.',fixed=T)        
        varstemp[1]=as.numeric(tempScan[[1]][2])#scan --good
        varstemp[2]=hit1[i,7]  ##tempScan[[1]][4] #charge --good
        varstemp[3]=xyCoord[[i]][1,1] #exp mass  ##Needs to be fixed?
        #print(varstemp[3])
        varstemp[4]=hit1[i,5] # ions num  --good
        varstemp[5]=hit1[i,6] # ions den  --good
        varstemp[6]=hit1[i,4] # sp  ##this is the sp for the first hit, the other code uses the sp from the 10th hit-- dont know why.
        varstemp[7]=hit1[i,3] # xcorr for hit 1  --good
        varstemp[8]=hit2[i,3] # xcorr for hit 2  --good
        
        ############################## xcorrPrime calculation start
        
        varstemp[9]=10*log(varstemp[7])/(log(2*varstemp[2]*(hit1[i,9]-1))) #xcorrprime  ##AA in paper AA-1 in code?? --good
        
        ############################ xcorrPrime calculation end
        varstemp[10]=hit2[i,2]*10 #delta cn hit2  --good
        varstemp[11]=hit1[i,1] # mhmass(1) --good
        #print(varstemp[11])
        varstemp[12]=10*hit1[i,5]/hit1[i,6] #ionsratio    ## 100* in paper 10* in code?? --good
        varstemp[13]=hit1[i,9] #AA (number of amino acids) --good
        
        ############################ delta mass 1 and delta mass 2 calc start
        
        varstemp[14]=(1e6)*abs((xyCoord[[i]][1,1]-hit1[i,1])/xyCoord[[i]][1,1]) #dmass1      --good!!
        #print(varstemp[14])
        #bob=(1e6)*abs((xyCoord[[i]][1,1]-hit1[i,1])/hit1[i,1])
        #print(bob)
        varstemp[15]=(1e6)*abs((xyCoord[[i]][1,1]-hit2[i,1])/xyCoord[[i]][1,1]) #dmass2      --good!!
        
        ########################## delta mass 1 and delta mass 2 calc end
        
        #varstemp[15]=(hit1[i,3]-hit2[i,3])/hit1[i,3] # dc2 or delta c2
        
        ######################### start calculate number of K and R in sequence   #--good
        numK=length(lengthTemp[[1]])-sum(as.numeric(is.na(pmatch(lengthTemp[[1]],"K",duplicates.ok=TRUE))))
        numR=length(lengthTemp[[1]])-sum(as.numeric(is.na(pmatch(lengthTemp[[1]],"R",duplicates.ok=TRUE))))
              ## there is a cleaner way to do this using logic in the phosSites calc
        varstemp[16]=numK+numR-1 #KR  --good
        #print(varstemp[16])
              ## note that in the matlab code we subract 1 from this value; not sure if that is needed 
              ## here, but i think it is to account for the theoretical terminal K/R
        
        ######################## end calculate number of K and R in sequence
        
        ##################### ps pt py and phos sites calc start        
        
        varstemp[17]=sum(phosIndex) # PhosSites  --good but i dont think it actually is a variable in model
        numPS=0
        numPT=0
        numPY=0
        if (varstemp[17]!=0){
          for (p in 1:(varstemp[13]+sum(phosIndex))){
          #print(p)
              if (phosIndex[p]==1){
                  site=posTempSep[[1]][p-1]
                  #print(site)
                  if(site=="S"){
                     numPS=numPS+1
                     #print(site)
                  }
                  if(site=="T"){
                     numPT=numPT+1
                     #print(site)
                  }
                  if(site=="Y"){
                     numPY=numPY+1
                     #print(site)
                  }                
              }        
          }
        }
        varstemp[18]=numPS  #--good
        varstemp[19]=numPT  #--good
        varstemp[20]=numPY  #--good(assumed as previous 2 were good)
        
        #################### ps pt py and phos sites calc end
        
        #################### start number mean number fancymean median calc  ## see noise.m
        
        no=length(xyCoord[[i]][,2])-1 #subtract the peak for exp mass
        sorted=sort(xyCoord[[i]][2:leng(xyCoord[[i]]),2],T)#sort peaks note:F=ascending order
        #print(sorted)
        normalized=xyCoord[[i]][2:leng(xyCoord[[i]]),2]/max(xyCoord[[i]][2:leng(xyCoord[[i]]),2]) #--need to eliminate 1st entry
        xbarfull=mean(xyCoord[[i]][2:leng(xyCoord[[i]]),2])
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
        
        # WORKS! #print((isolate(xyCoord[[i]][2:length(xyCoord[[i]][,1]),])))
        
        loners=isolate(xyCoord[[i]])  #-good
        top=findTop(varstemp[13],loners[[1]]) #-good
                       
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
          
          ##temp fixes from fullvars.m
        if(mods[26]==80){mods[26]=79.966332}
        if(mods[3]==132.027){mods[3]=132.066069}
        if(mods[4]==146.043){mods[4]=146.081719}  
          ##end temp fixes
        
        pep=pepcal(pseq(posTemp),varstemp[2],mods) # --good
        print(pep) 
        MminusH=c(MHH3PO4,MH2H3PO4)
        #print(varstemp[17])
        #print(loners[[1]])
        #print(top)
        #print(pep)
        #print(MminusH)
        #print(varstemp[2])
        #print(varstemp[17])
        us=uscore(top,pep,MminusH,varstemp[2],varstemp[17])
        #print(us)
        #print(length(us))
        
        #################### noa nda nsa toa tda tsa etc... calc end
        
        vars[[i]]=cbind(as.matrix(us),t(as.matrix(varstemp)))
        }
        close(out1)
        close(dta1)

}
#v=validate(outfiles,vars)
#print('here')

return(vars)

}



            
            
            
            
            
            
            
            
            