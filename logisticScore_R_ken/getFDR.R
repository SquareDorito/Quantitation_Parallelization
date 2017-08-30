#this method generates the percentage of false positives at each 
#score level
getFDR <- function(v,score){
#xy=matrix(0,101,2)
#for(i in 1:101){
#	xy[i,1]=(i-1)/100
#	posAtLevel=score>=((i-1)/100)
##   print(posAtLevel)
##	print(as.numeric(v&posAtLevel))
#	tp=sum(as.numeric(v&posAtLevel))
#	fp=2*sum(as.numeric(!v&posAtLevel))
##	print(v)
##	print(posAtLevel)
##	print(tp)
##	print(fp)
#	xy[i,2]=fp/length(v)#fp/(fp+tp)
#}
##xy[,2]=rev(xy[,2])
#return(xy) 	
o=rev(order(score))
fp=2*cumsum(v[o])
tp=cumsum(!v[o])
xy=cbind(score[o],v[o],fp,tp,fp/(1:length(fp))*100,tp/length(tp)*100)
colnames(xy)<-c('Score','RevDB','fp','tp','fpr','%tp')
#print(score[o])
#print(o)
#print(score)
return(xy)
}