combineDtaOut<-function(dtaDir,parDir){
setwd(dtaDir)
dtafiles <- list.files(pattern=".txt")
#print(dtafiles)
totalfiles=length(dtafiles)
for(i in 1:totalfiles){
  dta1 <- file(description=dtafiles[i], open="rt")
  dta2 <- readLines(dta1)
  dtafiles[i]=strsplit(dtafiles[i],'.',fixed=T)[[1]]
  if(i%%4==0){
    write(dta2,file=paste("E:/Temp/RLogisticScore/DtaOut/",dtafiles[i],".dta",sep=""))
    close(dta1)}
  else if(i%%4==1){
    write(dta2,file=paste("E:/Temp/RLogisticScore/DtaOut1/",dtafiles[i],".dta",sep=""))
    close(dta1)}
  else if(i%%4==2){
    write(dta2,file=paste("E:/Temp/RLogisticScore/DtaOut2/",dtafiles[i],".dta",sep=""))
    close(dta1)}
  else if(i%%4==3){
    write(dta2,file=paste("E:/Temp/RLogisticScore/DtaOut3/",dtafiles[i],".dta",sep=""))
    close(dta1)}
}

setwd(parDir)
parfiles <- list.files(pattern=".txt")
#print(parfiles)
totalfiles=length(parfiles)
for(i in 1:totalfiles){
  par1 <- file(description=parfiles[i], open="rt",encoding='utf-16')
  par2 <- readLines(par1)
  if(i%%4==0){
    write(par2,file=paste("E:/Temp/RLogisticScore/DtaOut/",parfiles[i],sep=""))
    close(par1)}
  else if(i%%4==1){
    write(par2,file=paste("E:/Temp/RLogisticScore/DtaOut1/",parfiles[i],sep=""))
    close(par1)}
  else if(i%%4==2){
    write(par2,file=paste("E:/Temp/RLogisticScore/DtaOut2/",parfiles[i],sep=""))
    close(par1)}
  else if(i%%4==3){
    write(par2,file=paste("E:/Temp/RLogisticScore/DtaOut3/",parfiles[i],sep=""))
    close(par1)}
}

}
