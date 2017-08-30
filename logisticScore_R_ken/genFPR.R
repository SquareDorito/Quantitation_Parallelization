genFPR<-function(valDir){
	FPRdat=NULL
	setwd(valDir)
	valFiles<-list.files()
	for(i in 1:length(valFiles)){
		val1<-file(description=valFiles[i],open='rt')
		val2<-readLines(val1)
		val2=strsplit(val2,"\t",fixed=T)[[1]]
		temp=cbind(val2[1],val2[2])
		FPRdat=rbind(FPRdat,temp)
		close(val1)
	}
	rRank=rev(order(FPRdat[,1]))
	FPRdat[rRank,2]=2*cumsum(as.numeric(FPRdat[rRank,2]))
	FPRdat=cbind(FPRdat[rRank,],as.numeric(FPRdat[rRank,2])/length(valFiles))
	colnames(FPRdat)<-c('Score','2xCumSum(revDB)','FPR')
	write.csv(FPRdat,file="C:/Temp/FPR.csv",quote=F,row.names = F)	
}