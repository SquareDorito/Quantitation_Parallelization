filterFile<-function(AAseq,Xcorr,charge){# does what findAndFilterOut.pl does on a file by file basis

  if((charge==1)&&(Xcorr<1.5)){return(0)}
  if((charge==2)&&(Xcorr<2)){return(0)}
  if((charge==3)&&(Xcorr<2.5)){return(0)}
  splitSeq=strsplit(AAseq,"",T)
  #print(splitSeq[[1]])
  if(min(as.numeric(is.na(pmatch(c("B","J","O","U","X","Z"),splitSeq[[1]]))))==FALSE){return(0)}
  return(1)
}