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
  retVars=predictMain(1,"C:/Users/knoh1/Documents/Quantitation_Parallelization/logisticScore_R_ken/DtaOut")
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
    setwd("C:/Users/knoh1/Documents/Quantitation_Parallelization/logisticScore_R_ken/validations")
    write(paste(valDat[i,1],'\t',valDat[i,2],sep=""),file=Fnames[i])
  }
}

source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//getDMods.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//getSMods.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//pseq.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//outmodsFile.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//getmass.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//getavgmass.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//vectorFind.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//getNumAfter.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//isolate.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//pepcalc.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//findTop.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//uscore.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//filterFile.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//validate.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//makeValidArray.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//parse_dt_out2.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//matComp.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//predictMain.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//PredictGui.r')
source('C://Users//knoh1//Documents//Quantitation_Parallelization//logisticScore_R_ken//logreg2.R')
dtaDir='E://Temp//dta'
parDir='E://Temp//LogisticScoreInput'
outDir='E://Temp//out'
#combineDtaOut(dtaDir,parDir)    ##need to clear folder \\proteome\Data\LogisticTempFiles\DtaOut and ...\Validations first!
#q(save="no")
commandMain('C:\\Users\\knoh1\\Documents\\Quantitation_Parallelization\\logisticScore_R_ken\\reduced.RData',dtaDir,outDir)
