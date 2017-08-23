pseq<-function(pep){ #pseq.pl

pep2=strsplit(pep,'',fixed=T)
newPep=NULL
i=1
for(j in 1:length(pep2[[1]])){
  if(pep2[[1]][j]=='*'){
    if(newPep[i-1]==17){newPep[i-1]=16}
    else if(newPep[i-1]==19){newPep[i-1]=18}
    else if(newPep[i-1]==23){newPep[i-1]=22}
  }
  else if(pep2[[1]][j]=='A'){newPep[i]=1}
  else if(pep2[[1]][j]=='C'){newPep[i]=2}
  else if(pep2[[1]][j]=='D'){newPep[i]=3}
  else if(pep2[[1]][j]=='E'){newPep[i]=4}
  else if(pep2[[1]][j]=='F'){newPep[i]=5}
  else if(pep2[[1]][j]=='G'){newPep[i]=6}
  else if(pep2[[1]][j]=='H'){newPep[i]=7}
  else if(pep2[[1]][j]=='I'){newPep[i]=8}
  else if(pep2[[1]][j]=='K'){newPep[i]=9}
  else if(pep2[[1]][j]=='L'){newPep[i]=10}
  else if(pep2[[1]][j]=='M'){newPep[i]=11}
  else if(pep2[[1]][j]=='N'){newPep[i]=12}
  else if(pep2[[1]][j]=='P'){newPep[i]=13}
  else if(pep2[[1]][j]=='Q'){newPep[i]=14}
  else if(pep2[[1]][j]=='R'){newPep[i]=15}
  else if(pep2[[1]][j]=='S'){newPep[i]=17}
  else if(pep2[[1]][j]=='T'){newPep[i]=19}
  else if(pep2[[1]][j]=='V'){newPep[i]=20}
  else if(pep2[[1]][j]=='W'){newPep[i]=21}
  else if(pep2[[1]][j]=='Y'){newPep[i]=23}
  #print(newPep)
  if(pep2[[1]][j]!='*'){i=i+1}
}
return(newPep)
}