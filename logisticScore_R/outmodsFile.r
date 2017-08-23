outmodsFile<-function(str1){

mods=array(0,26)

mods[1]=getNumAfter(str1,"A=")
mods[2]=getNumAfter(str1,"C=")
mods[3]=getNumAfter(str1,"D=")
mods[4]=getNumAfter(str1,"E=")
mods[5]=getNumAfter(str1,"F=")
mods[6]=getNumAfter(str1,"G=")
mods[7]=getNumAfter(str1,"H=")
mods[8]=getNumAfter(str1,"I=")
mods[9]=getNumAfter(str1,"K=")
mods[10]=getNumAfter(str1,"L=")
mods[11]=getNumAfter(str1,"M=")
mods[12]=getNumAfter(str1,"M=")
mods[13]=getNumAfter(str1,"P=")
mods[14]=getNumAfter(str1,"Q=")
mods[15]=getNumAfter(str1,"R=")
mods[16]=getNumAfter(str1,"S*=")
mods[17]=getNumAfter(str1,"S=")
mods[18]=getNumAfter(str1,"T*=")
mods[19]=getNumAfter(str1,"T=")
mods[20]=getNumAfter(str1,"V=")
mods[21]=getNumAfter(str1,"W=")
mods[22]=getNumAfter(str1,"Y*=")
mods[23]=getNumAfter(str1,"Y=")
mods[24]=getNumAfter(str1,"Nterm-pep=") ## see getNumAfter.r
mods[25]=getNumAfter(str1,"Cterm-pep=") ## see getNumAfter.r
mods[26]=getNumAfter(str1,"STY*")
#mods[23]=getNumAfter(str1,"")##need 3 more c-terminus, n-terminus, and sty*

return(mods)

}