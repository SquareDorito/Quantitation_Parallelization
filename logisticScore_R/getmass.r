#this method returns the monoisotopic mass of amino acids.  
getmass<-function(){  #getmass.m
mass=NULL
mass[1]=71.037114#A
mass[2]=103.009184#C
mass[3]=115.026944#D
mass[4]=129.042594#E
mass[5]=147.068414#F
mass[6]=57.021464#G
mass[7]=137.058912#H
mass[8]=113.084064#I
mass[9]=128.094963#K
mass[10]=113.084064#L
mass[11]=131.040484#M
mass[12]=114.042928#N
mass[13]=97.052764#P
mass[14]=128.058578#Q
mass[15]=156.101111#R
#mass[16]=87.032029#S*
mass[16]=87.032029#S
#mass[18]=101.047679#T*
mass[17]=101.047679#T
mass[18]=99.068414#V
mass[19]=186.079313#W
#mass[22]=163.063329#Y*
mass[20]=163.063329#Y
mass[21]=0#C-term
mass[22]=0#N-term
mass[23]=113.084064#X -- L or I
mass[24]=114.534936#B--avg of N & D
mass[25]=128.550586#Z avg of Q & E
names(mass)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","+Cterm-pep","+Nterm-pep","X","B","Z")
return(mass)
}