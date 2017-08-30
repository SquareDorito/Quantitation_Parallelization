#this method stores the static modifications in a convenient way
getSMods<-function(ModsStr){
	if(ModsStr==""||ModsStr==" "){return(NULL)}
	splitmods=strsplit(ModsStr,split=" ",fixed=T)
	splitmods=splitmods[[1]]
	numMods=length(splitmods)
	sMods=matrix(0,numMods,1)
	lab=NULL
	sMods=length(splitmods)
	splitmods=strsplit(splitmods,split="=",fixed=T)
#	namesIndex=2*(1:numMods)-1
#	modsIndex=2*(1:numMods)
	for(i in 1:numMods){
		sMods[i]=as.numeric(splitmods[[i]][2])
		lab=c(lab,splitmods[[i]][1])
	}
	names(sMods)<-lab
	return(sMods)
	
}

