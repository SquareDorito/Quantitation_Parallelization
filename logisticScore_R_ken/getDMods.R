#this method stores the dynamic modifications in a convenient way
getDMods<-function(ModsStr){

	if(ModsStr==""||ModsStr==" "){return(NULL)}
	splitmods=strsplit(ModsStr,split="|",fixed=T)
	splitmods=splitmods[[1]]
	numMods=length(splitmods)
	dMods=matrix(0,numMods,1)
	lab=NULL
	dMods=length(splitmods)
	splitmods=strsplit(splitmods,split=" ",fixed=T)
	for(i in 1:numMods){
		dMods[i]=as.numeric(splitmods[[i]][2])
		lab=c(lab,splitmods[[i]][1])
	}
	names(dMods)<-lab
	return(dMods)	

}
