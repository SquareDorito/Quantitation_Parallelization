keybd<-function(key){
  #cat(key)
  if(key=="0"||key=="1"){return(as.numeric(key))}
  else{return(noquote("Incorrect Key"))}
}

validate<-function(outFiles,vars){
print(noquote("Please validate the following files with 1 or 0"))
v=matrix(-1,length(outFiles),1)
index=1
while(index<=length(outFiles)){
    if(length(vars[[index]])==2){index=index+1}
    else{
      val=getGraphicsEvent(prompt=paste(outFiles[index],": "),onKeybd=keybd)
      print(val)
      if(val==1||val==0){
        v[index]=val
        index=index+1
      }
    }
    

}

return(v)
}