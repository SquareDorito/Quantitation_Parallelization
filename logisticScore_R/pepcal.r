pepcal<-function(pep,charge,mods){   # correct?

mmass=getmass()
amass=getavgmass()

for(k in 1:23){
  if(mods[k]!=0){
    mmass[k]=mods[k]
    amass[k]=mods[k]
  }
}
phosup=array(0,23)
phosup[16]=1
phosup[18]=1
phosup[22]=1
phosup=phosup*mods[26]
mmass=mmass+phosup
amass=amass+phosup

Nterm=mods[24]
Cterm=mods[25]

hplus=1.007825 - 0.0005485798959;
oxygen=18.04474;

z=c(0,0,0)

theo=NULL
term=length(pep)
for(n in 1:charge){
if(n<=2){m=mmass}
else{m=amass}
#print(m)   #########ERROR because "term" is not the "seq" in pepcalc.m
b=array(0,term)
y=array(0,term)
b[1]=((m[pep[1]]+Nterm)/n)+hplus
y[term]=(m[pep[term]]+Cterm+oxygen)/n+hplus
for(k in 2:term){
  b[k]=b[k-1]+(m[pep[k]]/n)
  y[term-k+1]=y[term-k+2]+(m[pep[term-k+1]]/n)
  #print(y)
}
b[term]=b[term]+(Cterm+oxygen)/n
y[1]=y[1]+(Nterm/n)
theo=rbind(theo,cbind(b,pep,y),z)
}

return(theo)
}