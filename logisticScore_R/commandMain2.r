commandMain2<-function(modelDir,dtaDir,outDir){
#combineDtaOut(dtaDir,outDir)
load(modelDir)
retVars=predictMain(1,"E:/Temp/RLogisticScore/DtaOut2")
vars=retVars[[1]]
Fnames=retVars[[2]]
vars=as.data.frame(vars)
names(vars)<-c("scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
val=predict(coefficients(reduced),vars)
#print(vars)
#print(coefficients(reduced ))
#print(val)
valDat=cbind(val,vars[,38])
valDat=as.data.frame(valDat)
names(valDat)<-c("v","rev_db")
#rRank=rank(valDat[,2])
#valDat[rRank,2]=2*cumsum(valDat[rRank,2])
#write.table(c(valDat[,1],valDat[,2],valDat[,2]/length(valDat[,2])),file="E:/Temp/Validations/FPR.txt",quote=F,sep="\t",row.names=F,col.names=F)
#print(valDat)
#print(Fnames)
for(i in 1:length(Fnames)){
#print(paste(valDat[i,1],"\t",valDat[i,2],sep=""))
	setwd("E:/Temp/Validations/")
	write(paste(valDat[i,1],'\t',valDat[i,2],sep=""),file=Fnames[i])
}
}
