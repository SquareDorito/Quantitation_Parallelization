combineDtaOut<-function(dtaDir,parDir){
  setwd(dtaDir)
  dtafiles <- list.files(pattern=".txt")
  #print(dtafiles)
  totalfiles=length(dtafiles)
  for(i in 1:totalfiles){
    dta1 <- file(description=dtafiles[i], open="rt")
    dta2 <- readLines(dta1)
    dtafiles[i]=strsplit(dtafiles[i],'.',fixed=T)[[1]]
    write(dta2,file=paste("E:/Temp/RLogisticScore/DtaOut/",dtafiles[i],".dta",sep=""))
    close(dta1)
  }
  setwd(parDir)
  parfiles <- list.files(pattern=".txt")
  #print(parfiles)
  totalfiles=length(parfiles)
  for(i in 1:totalfiles){
    par1 <- file(description=parfiles[i], open="rt",encoding='utf-16')
    par2 <- readLines(par1)
    write(par2,file=paste("E:/Temp/RLogisticScore/DtaOut/",parfiles[i],sep=""))
    close(par1)
  }
}

commandMain<-function(modelDir,dtaDir,outDir){
  load(modelDir)
  retVars=predictMain(1,"E:/Temp/RLogisticScore/DtaOut")
  vars=retVars[[1]]
  Fnames=retVars[[2]]
  vars=as.data.frame(vars)
  names(vars)<-c("scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
  val=predict(coefficients(reduced),vars)
  valDat=cbind(val,vars[,38])
  valDat=as.data.frame(valDat)
  names(valDat)<-c("v","rev_db")
  for(i in 1:length(Fnames)){
    #print(paste(valDat[i,1],"\t",valDat[i,2],sep=""))
    setwd("E:/Temp/Validations/")
    write(paste(valDat[i,1],'\t',valDat[i,2],sep=""),file=Fnames[i])
  }
}

source('E:/Temp/RLogisticScore/getDMods.r')
source('E:/Temp/RLogisticScore/getSMods.r')
source('E:/Temp/RLogisticScore/pseq.r')
source('E:/Temp/RLogisticScore/outmodsFile.r')
source('E:/Temp/RLogisticScore/getmass.r')
source('E:/Temp/RLogisticScore/getavgmass.r')
source('E:/Temp/RLogisticScore/vectorFind.r')
source('E:/Temp/RLogisticScore/getNumAfter.r')
source('E:/Temp/RLogisticScore/isolate.r')
source('E:/Temp/RLogisticScore/pepcalc.r')
source('E:/Temp/RLogisticScore/findTop.r')
source('E:/Temp/RLogisticScore/uscore.r')
source('E:/Temp/RLogisticScore/filterFile.r')
source('E:/Temp/RLogisticScore/validate.r')
source('E:/Temp/RLogisticScore/makeValidArray.r')
source('E:/Temp/RLogisticScore/parse_dt_out2.r')
source('E:/Temp/RLogisticScore/matComp.r')
source('E:/Temp/RLogisticScore/predictMain.r')
source('E:/Temp/RLogisticScore/PredictGui.r')
source('E:/Temp/RLogisticScore/logreg2.R')
dtaDir='E:/Temp/dta'
parDir='E:/Temp/LogisticScoreInput'
outDir='E:/Temp/out'
combineDtaOut(dtaDir,parDir)    ##need to clear folder \\proteome\Data\LogisticTempFiles\DtaOut and ...\Validations first!
#q(save="no")
commandMain('E:/Temp/RLogisticScore/Reduced.Rdata',dtaDir,outDir)
