#this is an internal calculation used ultimately in the uscore
#method.  Basically we pick out certain important peaks and 
#normalize the intensity                          
isolate <- function(dta,e){
  #WORKS-- DO NOT EDIT   
## will be the equivalent of the file isolate.m 
  m=max(dta[,2])
  dta[,2]=dta[,2]/m
  loners=NULL
  len=length(dta[,1])
  e=matrix(0,1,len)
  #print(dta)
  for(i in 1:len){
    flag=0
    intensity=dta[i,2]   
    interval=vectorFind((dta[,1]<dta[i,1]+5)&(dta[,1]>dta[i,1]-5))
    
    this=vectorFind(interval==i)
    if(this>1){
      leftpeak=interval[this-1]
      if((dta[i,1]-dta[leftpeak,1]<=1.2)&&(dta[i,1]-dta[leftpeak,1]>=.8)&&(dta[leftpeak,2]>=(1/2)*intensity)){
        flag=1
        e[i]=1
      }  
    }
    
    # the following lines followed by "#" correspond to the matlab command: interval[this]=[]
    index=1:length(interval)#
    indexRm=array(0,length(interval))#
    indexRm[this]=this#
    interval=interval[!indexRm&index]#
    interval10=vectorFind((dta[,1]<dta[i,1]+10)&(dta[,1]>dta[i,1]-10))
    this=vectorFind(interval10==i)
    # the following lines followed by "#" correspond to the matlab command: interval10[this]=NULL
    index=1:length(interval10)#
    indexRm=array(0,length(interval10))#
    indexRm[this]=this#
    interval10=interval10[!indexRm&index]#
    if(length(interval)>0){
      #range10=dta[interval,] #note this is commented out becasue the next line does not seem to work with sum(range10[,2])
      intensum=sum(dta[interval,2])
      if(dta[i,2]<1.5*(intensum/length(interval))){
        flag=1
        e[i]=2
      }
      for(j in 1:length(interval)){
        if(dta[interval[j],2]>intensity){
          flag=1
          e[i]=3
          break
        }
        else{
          spreaddist=abs(interval[j]-i)
          if(spreaddist<=2){
            dist=abs(dta[interval[j],1]-dta[i,1])
            if((dist<=spreaddist+.2)&&(dist>=spreaddist-.2)&&(dta[interval[j],2]>(3/4)*dta[i,2])){
              flag=1
              e[i]=4
              break
            }
            
          }
          
        }
      
      }
      
    }
  if(flag==0)
    loners=rbind(loners,dta[i,])  
    
  }
retvars=NULL  
retvars[[1]]=loners 
retvars[[2]]=e

return(retvars) 
}    
