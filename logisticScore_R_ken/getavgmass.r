#this method returns the average mass of amino acids.  
getavgmass<-function(){  #getavgmass.m
mass=NULL
mass[1]=71.0788#A
mass[2]=103.1388#C
mass[3]=115.0886#D
mass[4]=129.1155#E
mass[5]=147.1766#F
mass[6]=57.0519#G
mass[7]=137.1411#H
mass[8]=113.1594#I
mass[9]=128.1741#K
mass[10]=113.1594#L
mass[11]=131.1926#M
mass[12]=114.1038#N
mass[13]=97.1167#P
mass[14]=128.1307#Q
mass[15]=156.1875#R
#mass[16]=87.0782#S*
mass[16]=87.0782#S
#mass[18]=101.1051#T*
mass[17]=101.1051#T
mass[18]=99.1326#V
mass[19]=186.2132#W
#mass[22]=163.1760#Y*
mass[20]=163.1760#Y
mass[21]=0#C-term
mass[22]=0#N-term
mass[23]=113.1594#X -- L or I
mass[24]=114.5962#B--avg of N & D
mass[25]=128.6231#Z avg of Q & E
names(mass)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","+Cterm-pep","+Nterm-pep","X","B","Z")
return(mass)
}