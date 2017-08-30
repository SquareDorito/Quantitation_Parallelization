vectorFind <- function(x){
#print(x)
x=as.numeric(x)
## will be the equivalent of matlab "find"
z=NULL
count=1
for(i in 1:length(x)){
    if(x[i]==1){
      z[count]=i
      count=count+1
      }
}
return(z)

}
