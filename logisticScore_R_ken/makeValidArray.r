makeValidArray<-function(v,vars){

  newV=NULL
  count=1
  for(i in 1:length(vars)){
    if(length(vars[[i]])>2){
      newV[count]=vars[[i]]
      count=count+1
    }
  }
  return(newV)
}