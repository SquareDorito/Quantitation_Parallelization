#in this file we calculate the following variables:
#   1   Number of assigned peaks over number of peaks in top
#   2   Number of direct hits over number of peaks in top
#   3   Number of strong or direct hits over number of peaks in top
#   4   Total intensity of assigned peaks over total intensity of peaks in
#   top
#   5   Total intensity of direct hits over total intensity of peaks in top
#   6   Total intensity of strong or direct hits over total intensity of
#   peaks in top
#   7   Number of unassigned amino acids over peptide length
#   8   Number of amino acids not directly or strongly assigned over peptide
#   length
#   9   Number of amino acids not directly assigned over peptide length
#   10  Number of amino acids with more than one hit over peptide length
#   11  Number of amino acids with more than one strong or direct hit over peptide
#   length
#   12  Number of amino acids with more than one direct hit over peptide
#   length
leng<-function(vec){  ###need to correct Rbind and Cbind -- ensure "matrix" arguments
options(warn=-1)
a=max(ncol(vec),nrow(vec))
if(is.na(a)||is.null(ncol(vec))){return(length(vec))}
else{return(a)}
}

uscore<-function(top,pepca,MminusH,charge,phosSite,phosSymbol){
#print(top)
#print('-------')
#print(pepca)
#print('-------')
#print(MminusH)
#print('-------')
#print(charge)
#print('-------')
#print(phosSite)
topLength=leng(top[,1]) #length?     --good
pbreaks=vectorFind(pepca[,1]==0) #--good
aa=pbreaks[1]-1
#print((pbreaks))
for(k in 1:leng(pbreaks)){
  if(k==1){
    p=as.matrix(pepca[1:aa,1])
    #print(p)
    #print(pepca)
    p=rbind(p,as.matrix(pepca[1:aa,3]))
  }
  else{
  #print(pepca)
  #print('-------')
  #print(as.matrix(pepca))
  #print('-------')
  #print(pbreaks)
  #print('-------')
  #bob=(pbreaks[k-1]+1):(pbreaks[k]-1)
  #print(bob)
  #rint('-------')
    p=rbind(p,as.matrix(pepca[(pbreaks[k-1]+1):(pbreaks[k]-1),1]))
    p=rbind(p,as.matrix(pepca[(pbreaks[k-1]+1):(pbreaks[k]-1),3]))
  }
}
#print(pepca)
#print(leng(pepca))
p=rbind(p,as.matrix(pepca[leng(pepca)-1,1]),0)        ## length?   --PREVIOUS CODE IS CORRECT
################# GOOD #################
#print(p)

if(phosSite>0){
  p[leng(p)]=pepca[leng(pepca)-1,1]-(98/charge)    ##length?
}
#print(p)
pep=pepca[1:aa,2]
#print(pep) #good till here   
m=0
for(i in 2:(aa+1)){
  if((pep[i-1]==paste("S",phosSymbol,sep=""))||(pep[i-1]==paste("T",phosSymbol,sep=""))){
    m[i]=1
  }
  else{
    m[i]=m[i-1]
  }
}
m=m[2:(aa+1)]
#print(m)
####################
my=0
for(i in 2:(aa+1)){
	#print(pep[i-1])
  j=aa+3-i
  if((pep[j-1]==paste("S",phosSymbol,sep=""))||(pep[j-1]==paste("T",phosSymbol,sep=""))){
    my[j]=1
  }
  else{
	if(j==aa+1){
		my[j]=0
	}else{
		my[j]=my[j+1]
	}
  }
}
my=my[2:(aa+1)]
####################
###################################### above should be correct!!            ##
divvec=NULL
mvec=NULL
endsvec=NULL
blackout=array(1,aa)
blackout[aa]=0

for(n in 1:charge){
  divvec=rbind(divvec,n*matrix(1,2*aa,1))
  #mvec=rbind(mvec,as.matrix(m),matrix(m[leng(m):1]))
  mvec=rbind(mvec,as.matrix(m),as.matrix(my))
  #print(leng(blackout))
  endsvec=rbind(endsvec,matrix(blackout),matrix(blackout[leng(blackout):1]))
}
divvec=rbind(matrix(divvec),charge*matrix(1,2,1))
mvec=rbind(mvec,matrix(0,2,1))
if(phosSite>1){
  mvec[leng(mvec)]=1
}
#print(mvec) --GOOD
endsvec=rbind(matrix(endsvec),matrix(1,2,1))

wloss=18/divvec
aloss=17/divvec
mloss=98/divvec
                                 
#########LINE 105 uscore.m    GOOD TILL HERE. 
#print(p)#[,1])#-wloss)
#print(endsvec)

#print(p)
#print('-------')
p=cbind(p,(p[,1]-wloss)*endsvec)
p=cbind(p,(p[,1]-aloss)*endsvec)
p=cbind(p,(p[,1]-2*wloss)*endsvec)
p=cbind(p,(p[,1]-2*aloss)*endsvec)
p=cbind(p,(p[,1]-mloss)*mvec*endsvec)

p[leng(p)-1,1]=0  ## length?
#good.

inds=matrix(0,leng(p[,1]),6)
ints=matrix(0,leng(p[,1]),6)

top=cbind(top,matrix(0,leng(top[,1]),1))
#print(top)
ups=matrix(1:3,3,1)
for(b in 1:leng(p[,1])){
  for(d in 1:6){
    f=vectorFind((top[,1]<p[b,d]+.5)&(top[,1]>p[b,d]-.5))
    u=ups
    if(b==(leng(p[,1])-1)){u=u+1}
    if(length(f)>0){ #line 130
      inds[b,d]=1
      pks=rbind(top[f,])
      ints[b,d]=max(pks[,2])
      mz=vectorFind((pks[,2])==ints[b,d])
      for(za in 1:length(mz)){
        spot=vectorFind(top[,1]==pks[mz[za]])
        if(d==1){top[spot,3]=u[3]}
        else if(d==2||d==3){top[spot,3]=max(top[spot,3],u[2])}
        else{top[spot,3]=max(top[spot,3],u[1])}
      }
    }
  }
}

##line 150
maxindloss=max(inds[length(p[,1])-1,2:3])     #presumed correct
maxintloss=max(ints[length(p[,1])-1,2:3])     #presumed correct.
inds[length(p[,1])-1,]=c(maxindloss,inds[length(p[,1])-1,4:5],0,0,0)
ints[length(p[,1])-1,]=c(maxintloss,ints[length(p[,1])-1,4:5],0,0,0)
#good
hits=matrix(0,1,4)
intens=matrix(0,1,4)
for(kk in 1:length(top[,1])){
  hits[top[kk,3]+1]=hits[top[kk,3]+1]+1
  intens[top[kk,3]+1]=intens[top[kk,3]+1]+top[kk,2]
}
#print(hits)
#print(intens)#############################################################################################################

if(sum(hits[2:4])==0){
  nda=0
  ndi=0
  nsa=0
  nsi=0
}
else{
  nda=hits[4]/sum(hits[2:4])
  ndi=intens[4]/sum(intens[2:4])
  nsa=sum(hits[3:4])/sum(hits[2:4])
  nsi=sum(intens[3:4])/sum(intens[2:4])
}

if(sum(hits[1:4])==0){
noa=0
tai=0
}
else{
noa=sum(hits[2:4])/sum(hits[1:4])
tai=sum(intens[2:4])/sum(intens[1:4])
}
##line 186
countvars=cbind(noa,nda,nsa)
intenvars=cbind(tai,ndi,nsi)
#good.
indsum=NULL
indsum=cbind(indsum,(apply(inds,1,sum)))##line 190--correct?
indsum=cbind(indsum,(apply(inds[,1:3],1,sum)))
indsum=cbind(indsum,inds[,1])

intsum=NULL
intsum=cbind(intsum,(apply(ints,1,sum)))
intsum=cbind(intsum,(apply(ints[,1:3],1,sum)))
intsum=cbind(intsum,ints[,1])
#good.
##Line 202         --resume testing here.
Bdcoverage=NULL
Ydcoverage=NULL
fullcoverage=NULL
Btcoverage=NULL
Ytcoverage=NULL
zcount=NULL
tcount=NULL
#print(indsum)
#print(intsum)
for(u in 1:3){
  Bd=matrix(0,aa,1)
  Yd=matrix(0,aa,1)
  Bt=matrix(0,aa,1)
  Yt=matrix(0,aa,1) 
  for(v in 1:charge){
    addBd=indsum[((v-1)*2*aa+1):((v-1)*2*aa+aa),u]               ##error in this line
    addYd=indsum[((v-1)*2*aa+aa+1):(v*2*aa),u]
    #print('here') 
    Bd=cbind(Bd,matrix(addBd))
    Yd=cbind(Yd,matrix(addYd))
    addBt=intsum[((v-1)*2*aa+1):((v-1)*2*aa+aa),u]
    addYt=intsum[((v-1)*2*aa+aa+1):(v*2*aa),u]
    Bt=cbind(Bt,addBt)
    Yt=cbind(Yt,addYt)
  }                        ##matlab sum(arg)-->apply(arg,2,sum) 
  #print('here')
  Bdcoverage=cbind(Bdcoverage,apply(Bd,1,sum))# correct?  line 219
  Ydcoverage=cbind(Ydcoverage,apply(Yd,1,sum))
  fullcoverage=cbind(fullcoverage,Bdcoverage[,u]+Ydcoverage[,u])
  zcount=cbind(zcount,length(vectorFind(fullcoverage[,u]==0))/aa)
  tcount=cbind(tcount,length(vectorFind(fullcoverage[,u]>=2))/max(length(vectorFind(fullcoverage[,u]>=1)),1))
  Btcoverage=cbind(Btcoverage,apply(Bt,1,sum))
  Ytcoverage=cbind(Ytcoverage,apply(Yt,1,sum))
}
##line 230

coveragevars=cbind(zcount,tcount)   ##error in zcount!!  --corrected.
#print(coveragevars) #-- not the same as MATLAB
prolines=as.numeric(pepca[1:aa,2]=="P")
proMinus1=c(prolines[2:length(prolines)],0) #need b-ion @ this index
proPlus1=c(0,prolines[1:(length(prolines)-1)]) #need y-ion @ this index
#print(prolines)
#print(fullcoverage)
#if(sum(prolines)>0){
#  pro=matrix(0,sum(prolines),3)
#  proB=matrix(0,sum(prolines),3)
#  proYm1=matrix(0,sum(proPlus1),3)
#  proBm1=matrix(0,sum(proMinus1),3)
#  #proind=rbind(fullcoverage[prolines,])##Not sure how to ensure dimensions arent transposed-- this seems to work. 
#  #the following changed 6/3/08
#  proind=rbind(Ydcoverage[prolines,])
#  proBind=rbind(Bdcoverage[prolines,])
#  proYm1ind=rbind(Ydcoverage[proPlus1,])
#  proBm1ind=rbind(Bdcoverage[proMinus1,])
#  #print(proind)
#  #proint=rbind(Btcoverage[prolines,]+Ytcoverage[prolines,],matrix(0,1,3))
#  #the following changed 6/3/08
#  proint=rbind(Ytcoverage[prolines,],matrix(0,1,3))
#  proBint=rbind(Btcoverage[prolines,],matrix(0,1,3))
#  proYm1int=rbind(Ytcoverage[proPlus1,],matrix(0,1,3))
#  proBm1int=rbind(Btcoverage[proMinus1,],matrix(0,1,3))
#    #print(proint)-ok!
#  for(g in 1:3){
#    nonz=vectorFind(proind[,g]>0)
#    nonzB=vectorFind(proBind[,g]>0)
#    nonzYm1=vectorFind(proYm1ind[,g]>0)
#    nonzBm1=vectorFind(proBm1ind[,g]>0)
#    pro[nonz,g]=1
#    proB[nonzB,g]=1
#    proYm1[nonzYm1,g]=1
#    proBm1[nonzBm1,g]=1
#  }
#  pro=rbind(pro,matrix(0,1,3))
#  print(pro)
#  proB=rbind(proB,matrix(0,1,3))
#  proYm1=rbind(proYm1,matrix(0,1,3))
#  proBm1=rbind(proBm1,matrix(0,1,3))
#  profound=apply(pro,2,sum)/sum(prolines)
#  print(profound)
#  proBfound=apply(proB,2,sum)/sum(prolines)
#  proYm1found=apply(proYm1,2,sum)/sum(prolines)# use prolines to avoid by divide by 0 error-- would default to 0 anyway
#  proBm1found=apply(proBm1,2,sum)/sum(prolines)# use prolines to avoid by divide by 0 error-- would default to 0 anyway
#  propossible=leng(prolines)/aa
#  print(propossible)
#  proYm1possible=leng(proMinus1)/aa
#  proBm1possible=leng(proPlus1)/aa
#  avgproint=apply(proint,2,sum)/apply(rbind(apply(pro,2,sum),matrix(1,1,3)),2,max)    ##??
#  print(avgproint)
#  avgBproint=apply(proBint,2,sum)/apply(rbind(apply(proB,2,sum),matrix(1,1,3)),2,max)    ##??
#  avgYm1proint=apply(proYm1int,2,sum)/apply(rbind(apply(proYm1,2,sum),matrix(1,1,3)),2,max)    ##??
#  avgBm1proint=apply(proBm1int,2,sum)/apply(rbind(apply(proBm1,2,sum),matrix(1,1,3)),2,max)    ##??
#  #print(avgproint)
#  topavginten=apply(cbind(top[,2]),2,sum)/leng(top[,2])   ##?
#  normavgproint=avgproint/topavginten
#  print(normavgproint)
#  normavgBproint=avgBproint/topavginten
#  normavgYm1proint=avgYm1proint/topavginten
#  normavgBm1proint=avgBm1proint/topavginten
#  provars=cbind(t(as.matrix(profound[3])),t(as.matrix(proBfound[3])),t(as.matrix(proYm1found[3])),t(as.matrix(proBm1found[3])),propossible,proYm1possible,proBm1possible,t(as.matrix(avgproint[3])),t(as.matrix(avgBproint[3])),t(as.matrix(avgYm1proint[3])),t(as.matrix(avgBm1proint[3])),t(as.matrix(normavgproint[3])),t(as.matrix(normavgBproint[3])),t(as.matrix(normavgYm1proint[3])),t(as.matrix(normavgBm1proint[3])))
#  print(provars)
#}
#else{
#  provars=matrix(0,1,39)
#}
#line 260/300
#print(countvars)
#print(intenvars)
#print(coveragevars)
#print(provars)
#outs=cbind(countvars,intenvars,coveragevars,provars)
outs=cbind(countvars,intenvars,coveragevars)
return(outs)
}