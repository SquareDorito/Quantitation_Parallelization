#this method basically resorts the variables.  It was initially 
#writen to compare the R output with Matlab, but later it was just 
#easier to keep varibales in this order
matComp<-function(vars2,Fnames){
vars=NULL
Fnames2=NULL
count=1
for(j in 1:length(vars2)){
  if(length(vars2[[j]])>2){
    vars[[count]]=vars2[[j]]
    Fnames2[count]=Fnames[j]
    count=count+1
  }
}
newV=matrix(0,38,length(vars))
for(i in 1:length(vars)){
  newV[1,i]=vars[[i]][13]
  newV[2,i]=vars[[i]][14]
  newV[3,i]=vars[[i]][19]
  newV[4,i]=vars[[i]][20]
  newV[5,i]=vars[[i]][21]
  newV[6,i]=vars[[i]][22]
  newV[7,i]=vars[[i]][23]
  newV[8,i]=vars[[i]][16]
  newV[9,i]=vars[[i]][17]
  newV[10,i]=vars[[i]][24]
  newV[11,i]=vars[[i]][18]
  newV[12,i]=vars[[i]][25]
  newV[13,i]=vars[[i]][28]
  newV[14,i]=vars[[i]][29]
  newV[15,i]=vars[[i]][26]
  newV[16,i]=vars[[i]][27]
  newV[17,i]=vars[[i]][30]
  newV[18,i]=vars[[i]][31]
  newV[19,i]=vars[[i]][32]
  newV[20,i]=vars[[i]][33]
  newV[21,i]=vars[[i]][34]
  newV[22,i]=vars[[i]][35]
  newV[23,i]=vars[[i]][36]
  newV[24,i]=vars[[i]][37]
  newV[25,i]=vars[[i]][38]
  newV[26:37,i]=vars[[i]][1:12]
  newV[38,i]=vars[[i]][39]
  }
  retVars=NULL
  retVars[[1]]=newV
  retVars[[2]]=Fnames2
  return(retVars)
}



