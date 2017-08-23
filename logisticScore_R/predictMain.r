# this file calls the main calculating function (parse_dt_out2)
# and formats the output appropriately (using matcomp) for either 
# "validated" or "to be validated"
predictMain<-function(method,curdir,valFile){

 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/Anthony/final_combinations/Rparsing/"
 #curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/matlabandR/Uploaded 12-14/Files/AlphaCasein/acasein_FTMSLTMSMS_HRG2_081105_155721"
 #curdir = "//proteome/User Files/old users/Lisa/mcp5 data/ctrl_IgE_mcp5_5e8"
 
  if(method==1){ #validate after vars are computed
    setwd(curdir)
    
    retVars=parse_dt_out2()
    vars=retVars[[1]]
    Fnames=retVars[[2]]

    retV=matComp(vars,Fnames)
    varMat=retV[[1]]
    Fnames=retV[[2]]
    #print(retV)
    #print(Fnames)
    varMat=t(varMat[1:38,])
    #print(dir)

  }
  else if(method==2){   #validate from file   format: <1 or 0>"\t"<scan>"\t"<charge>
    setwd(curdir)
    #print(valFile)
    vars=parse_dt_out2()
    Fnames=vars[[2]]
    varMat=matComp(vars[[1]],vars[[2]])[[1]]
    varMat=rbind(varMat,-1)
    valid=read.delim(valFile,header=F)
    #print(valid)
#    print('here')
#    for(i in 1:length(valid[,1])){
#      print(varMat[1.])
#      print(valid[i,2])
	  varMat[1,]=varMat[1,]+varMat[2,]/10
	  order1=order(varMat[1,])
	  valid[,2]=valid[,2]+valid[,3]/10
	  print(valid)
	  order2=order(valid[,2])
	  varMat[39,order1]=valid[order2,1]
#      ind=vectorFind(varMat[1,]==valid[i,2])
#      #print(ind)
#      ind2=1
#      if(length(ind)>0){
#        ind2=vectorFind(varMat[2,ind]==valid[i,3])
#        #print(ind[ind2])
#        if(length(ind2)>0){
#          varMat[39,ind[ind2]]=valid[i,1]
#        }
#      }
    
#    }
  }
  #print(varMat)
  retVars=NULL
  retVars[[1]]=varMat
  retVars[[2]]=Fnames
  return(retVars)
}
