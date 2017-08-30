#this method finds every possible theoretical peak for a given 
#peptide (excluding possiblities of loss-- this is taken care of 
#later).  
pepcalc<-function(pep,charge,dmods,smods){

	dMods=getDMods(dmods)
	#print(dMods)
	sMods=getSMods(smods)
	amass=getavgmass()
	mmass=getmass()
	if(length(sMods>0)){
		for(i in 1:length(sMods)){
			amass[names(sMods)[i]]=sMods[names(sMods)[i]]
			mmass[names(sMods)[i]]=sMods[names(sMods)[i]]
		}
	}
	#print(amass)
	#print(mmass)
	splitPep=strsplit(pep,'',fixed=T)[[1]]
	pepMass=NULL
	apepMass=NULL
	seq=NULL
	for(p in 1:length(splitPep)){
		if(!is.na(mmass[splitPep[p]])){
			pepMass=c(pepMass,mmass[splitPep[p]])
			apepMass=c(apepMass,amass[splitPep[p]])
			seq=c(seq,splitPep[p])
		}
		else{
			pl=length(pepMass)
			#print(dMods[splitPep[p]])
#			print(class(dMods[splitPep[p]]))
#			print(splitPep[p])
#			print(dMods)
			pepMass[pl]=pepMass[pl]+dMods[splitPep[p]]
			apepMass[pl]=apepMass[pl]+dMods[splitPep[p]]
			seq[pl]=paste(seq[pl],splitPep[p],sep="")
		}
	}
	bIons_1=cumsum(pepMass)
	yIons_1=cumsum(rev(pepMass))
	last=length(bIons_1)
	bIons_1=bIons_1+mmass["+Nterm-pep"]+1.00727
	bIons_1[last]=bIons_1[last]+18.0105+mmass["+Cterm-pep"]
	yIons_1=yIons_1+18.0105+mmass["+Cterm-pep"]+1.00727
	yIons_1[last]=yIons_1[last]+mmass["+Nterm-pep"]	
	yIons_1=rev(yIons_1)
	Block=cbind(as.matrix(bIons_1),as.matrix(seq),as.matrix(yIons_1))
	#print(charge)
	if(charge>1){
		#pepMass=NULL
		#seq=NULL
		yIons=matrix(0,last,charge-1)
		bIons=matrix(0,last,charge-1)
		#for(p in 1:length(splitPep)){
		#if(!is.na(mmass[splitPep[p]])){
		#	pepMass=c(pepMass,mmass[splitPep[p]])
		#	seq=c(seq,splitPep[p])
		#}
		#else{
		#	pl=length(pepMass)
		#	pepMass[pl]=pepMass[pl]+dMods[splitPep[p]]
		#	seq[pl]=paste(seq[pl],splitPep[p],sep="")
		#}
		#}
	abIons_1=cumsum(apepMass)  #average
	ayIons_1=cumsum(rev(apepMass))
	last=length(abIons_1)
	abIons_1=abIons_1+amass["+Nterm-pep"]+1.00727
	abIons_1[last]=abIons_1[last]+18.0105+amass["+Cterm-pep"]
	ayIons_1=ayIons_1+18.0105+amass["+Cterm-pep"]+1.00727
	ayIons_1[last]=ayIons_1[last]+amass["+Nterm-pep"]	
	ayIons_1=rev(ayIons_1)
		for(i in 1:(charge-1)){
			#bIons[,i]=cumsum(pepMass)
			#yIons[,i]=cumsum(rev(pepMass))
			#bIons[,i]=bIons[,i]+mmass["+Nterm-pep"]
			#bIons[,i]=bIons[,i]/(i+1)+1.00727
				bIons[,i]=(abIons_1+1.00727*i)/(i+1)
			#bIons[last,i]=bIons[last,i]+18.0105+mmass["+Cterm-pep"]
			#yIons[,i]=yIons[,i]+18.0105+mmass["+Cterm-pep"]
			#yIons[last,i]=yIons[last,i]+mmass["+Nterm-pep"]
			#yIons[,i]=yIons[,i]/(i+1)+1.00727
				yIons[,i]=(ayIons_1+1.00727*i)/(i+1)
			#yIons[,i]=rev(yIons[,i])
			Block2=cbind(bIons[,i],as.matrix(seq),yIons[,i])
			#print(Block2)
			Block=rbind(Block,c(0,0,0),Block2)
		}
	}
	Block=rbind(Block,c(0,0,0))
	Block=data.frame(as.numeric(Block[,1]),Block[,2],as.numeric(Block[,3]))
	#print(Block)
	return(Block)
}
	