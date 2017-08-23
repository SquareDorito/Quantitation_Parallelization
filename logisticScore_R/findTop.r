findTop<-function(num,vec){### this method will find the "num" largest rows in vec based on col 2

#top=vec[1:num,]
#minTop=min(top[,2])
#minIndex=which.min(top[,2])
#for(i in (num+1):length(vec[,2])){
#  #print(length(vec))
#  #print(i)
#  if(vec[i,2]>minTop){
#    top[minIndex,]=vec[i,]
#    minTop=min(top[,2])
#    minIndex=which.min(top[,2])
#  }
#}

s1<-sort.list(vec[,2],decreasing = T)
temp=vec[s1,]
top=temp[1:num,]
#now check for duplicates in top[,2]--matlab defers to col 1 if duplicate, here we do the same thing
if(length(vec[,1])<num){
for(i in 1:num){
  if(temp[i,2]==temp[i+1,2]){  
    if(temp[i,1]<temp[i+1,1]){
      #print('here')
      temp2=temp[i,]
      top[i,]=temp[i+1,]
      if(i!=num){
        top[i+1,]=temp2
      }
    }  
  }
}
}

return(top)
}
