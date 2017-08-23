#this file contains all of the code required for the GUI 
#and makes all of the calls to the computational methods. 
#Moreover, the parser (written in perl) is called in this file.  

applyOldModel<-function(scriptLocation){
##find model you want to use
  tkmessageBox(message="Please find the .RData file with the model you want to use.")
  modelFile<-tkgetOpenFile(filetypes = "{{Model Files} {.RData}}")
  print(tclvalue(modelFile))
#    
#  modelFile2<-NULL
#  for(i in 1:length(modelFile)){
#    modelFile2=paste(modelFile2,modelFile[i])
#  }
#  load(sub("^ ","",modelFile2))
  load((tclvalue(modelFile)))
#  model<-as.character(modelFile2)
  model<-as.character(tclvalue(modelFile))
  model<-strsplit(model,"/",T)
  model<-strsplit(model[[1]][length(model[[1]])],".RData",T)
  model<-model[[1]][1]
##find dta/out files to apply to
  tkmessageBox(message="Please find the directory with the .dta and .out files you wish to apply the model to")
  curdir<-tkchooseDirectory()
  #######################################################################call perl script to pre-parse .out files
# setwd(scriptLocation )
# #system(paste("perl parseOut.pl",curdir))
# system(paste("perl parsepepdoc.pl '",tclvalue(curdir),"'",sep=""),ignore.stderr=T,wait=T)
 ##########################################################################################
#  curdir<-as.character(curdir)
#  curdir2=NULL
#  for(i in 1:length(curdir)){
#       curdir2=paste(curdir2,curdir[i])
#    } 
  vars=predictMain(1,sub("^ ","",tclvalue(curdir)))[[1]]
  #vars=matComp(vars)
  vars=as.data.frame(vars)
  colnames(vars)<-c("scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
  #print(vars)
  if(model=="sequest"){valid=predict(coefficients(sequest),vars)}
  else if(model=="sequestplus"){valid=predict(coefficients(sequestplus),vars)}
  else if(model=="full"){valid=predict(coefficients(full),vars)}
  else if(model=="reduced"){valid=predict(coefficients(reduced),vars)}
  else{print(noquote("cannot load model"))}
  validations=cbind(valid,vars[,1:2],vars[,38])
  colnames(validations)<-c('valid','scan','charge','reverse DB')
  write.csv(validations,file=paste(model,"validations.csv",sep=""))
  name=paste(model,"validations.csv")
  fpr=getFDR(as.numeric(as.logical(vars[,38])),validations[,1])
  name2=paste(model,"fpr.csv",sep="")
  write.csv(fpr,file=name2)
  tkmessageBox(message=noquote(paste("the validations are in the file\n",name,'\nthe FPR is in the file\n',name2,"\nin the directory\n",tclvalue(curdir),sep="")))
 
##write validations 

}
####################################################################################
makeNewModel<-function(scriptLocation){

  require(tcltk)
  tt<-tktoplevel()
  tkwm.title(tt,"Validate")
  done<-tclVar(0)
  Yes.but <- tkbutton(tt,text="     Yes     ",    command=function() tclvalue(done)<-1)
  No.but <- tkbutton(tt,text="     No      ",command=function() tclvalue(done)<-2)

  # Place the two buttons on the same row in their assigned window (tt).
  qfont<-tkfont.create(family="courier",size=16,weight="bold")
  tkgrid(tklabel(tt,text="  "),columnspan=4)
  tkgrid(tklabel(tt,text="Is Your Data Set Already Validated?",font=qfont),columnspan=2)
  tkgrid(tklabel(tt,text="  "))
  tkgrid(Yes.but,No.but)
  tkgrid(tklabel(tt,text="  "))
  tkwait.variable(done)
  doneVal <- as.integer(tclvalue(done))
  if(doneVal==1){
    tkdestroy(tt)
    tkmessageBox(message="Please find the directory with the .dta and .out files you wish to build model with")
    curdir<-tkchooseDirectory()
     #######################################################################call perl script to pre-parse .out files
 setwd(scriptLocation )
 #system(paste("perl parseOut.pl",curdir))
 system(paste("perl parsepepdoc.pl '",tclvalue(curdir),"'",sep=""),ignore.stderr=T,wait=T)
   # print("HERE")
##########################################################################################
    #print(curdir)
    #print(dir())
    tkmessageBox(message="Please find the tab-deliminated file with the validations")
    valFile<-tkgetOpenFile(filetypes = "{{Tab-Delim Files} {.txt .tab}}")
    #valFile<-as.character(valFile)
#    curdir<-as.character(curdir)
    curdir2=tclvalue(curdir)
    valFile2=tclvalue(valFile)
#    for(i in 1:length(valFile)){
#       valFile2=paste(valFile2,valFile[i])
#    }
#    for(i in 1:length(curdir)){
#       curdir2=paste(curdir2,curdir[i])
#    }
    #print(curdir2)
    #print(valFile2)

    varValMat=predictMain(2,sub("^ ","",curdir2),sub("^ ","",valFile2))[[1]]  ## NOW PLUG INTO REGRESSION METHOD
    #varValMat=(t(varValMat))
    varValMat=rbind(varValMat[39,],varValMat[1:38,])
    varValMat=as.data.frame(t(varValMat))
    #print(varValMat)
    names(varValMat)<-c("v","scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
    #print(varValMat)
    logreg2(sub("^ ","",curdir2),varValMat)
  }
  if(doneVal==2){     ##CONTINUE HERE.
    tkdestroy(tt)
    tkmessageBox(message="Please find the directory with the .dta and .out files you wish to build model with")
    curdir<-tkchooseDirectory()
    
    #######################################################################call perl script to pre-parse .out files
    setwd(scriptLocation )
    #system(paste("perl parseOut.pl",curdir))
 system(paste("perl parsepepdoc.pl '",tclvalue(curdir),"'",sep=""),ignore.stderr=T,wait=T)
    #print("HERE")
    ##########################################################################################
    #curdir<-as.character(curdir)
    curdir2=tclvalue(curdir)
#    for(i in 1:length(curdir)){
#       curdir2=paste(curdir2,curdir[i])
#    }
    Mat<-predictMain(1,sub("^ ","",curdir2))
    varValMat=(t(Mat[[1]]))
    l=length(Mat[[2]])
    #print(varValMat)
    #print(l)
    randSamp=sort(sample(1:l,max(75,round(l/10))))
    validateThese=cbind(-1,t(varValMat[,randSamp]))
    pfiles <- list.files(pattern=".txt")
    name=strsplit(pfiles[1],split=']',fixed=T)
    colnames(validateThese)<-c("v","scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
    #print(validateThese)
    write.csv(validateThese,file=paste(paste(name[[1]][1],"validate",sep="_"),".csv",sep=""))
    tkmessageBox(message=paste("Please validate",paste(paste(name[[1]][1],"validate",sep="_"),".csv",sep=""),"by replacing each -1 with either 1 or 0.  Press OK when done."))
    validateThese<-read.csv(paste(paste(name[[1]][1],"validate",sep="_"),".csv",sep=""))
    validateThese<-as.data.frame(validateThese[,2:40])
    

    #print(sub("^ ","",curdir2))
    #return(validateThese)
    varValMat=t(varValMat)
	colnames(varValMat)<-c("scan","charge","xcorr1","xcorr2","xcorrp","dc2","mhmass","ionsnum","ionsden","ionsratio","sp","aa","kr","phos","dmass1","dmass2","ps","pt","py","fancymean","number","median","mean","bscore","sumscore","noa","nda","nsa","toa","tda","tsa","percunass","percweakass","percnondirass","onehit","onestronghit","onedirecthit","rev_database")
    logreg2(sub("^ ","",curdir2),validateThese)
    #print(doneVal)
    #####the following should work.
    load(paste(sub("^ ","",curdir2),"/full.RData",sep=""))
    #print(coefficients(full))
    #print(varValMat)
    setwd(sub("^ ","",curdir2))
    valid=predict(coefficients(full),varValMat)
    #print(varValMat)
    validations=cbind(valid,varValMat[,1:2],varValMat[,38])
    colnames(validations)<-c('valid','scan','charge','reverse DB')
    write.csv(validations,file=paste(name[[1]][1],"FullValidations.csv",sep=""))
    name=paste(name[[1]][1],"FullValidations.csv",sep='')
    tkmessageBox(message=(paste("the validations are in the file\n",name,"\nin the directory\n",curdir2,'\nNote: Forward Database -> 0, Reverse Database -> 1',sep="")))
  }


}
####################################################################################
PredictGui<-function(scriptLocation){
  # Load the TclTk package
  require(tcltk)

  # Create a new toplevel window
  tt <- tktoplevel()

  # Give the window a title
  tkwm.title(tt,"Predictor")

  rb1 <- tkradiobutton(tt)
  rb2 <- tkradiobutton(tt)
  rbValue <- tclVar("oranges")
  tkconfigure(rb1,variable=rbValue,value="Make New Model")
  tkconfigure(rb2,variable=rbValue,value="Apply Old Model")
  fontHeading <- tkfont.create(family="courier",size=24,weight="bold")
  tkgrid(tklabel(tt,text="Predictor Main Function",font=fontHeading),columnspan=2)
  tkgrid(tklabel(tt,text="\t\t\t\t\t\t\t\t\t"),columnspan=2)
  tkgrid(tklabel(tt,text="What Would You Like To Do?"),columnspan=2)
  tkgrid(tklabel(tt,text="    "))
  tkgrid(tklabel(tt,text="Make New Model "),rb1)
  tkgrid(tklabel(tt,text="Apply Old Model "),rb2)
  tkgrid(tklabel(tt,text="    "))
  tkgrid(tklabel(tt,text="    "))
  OnOK <- function()
  {
      rbVal <- as.character(tclvalue(rbValue))

      tkdestroy(tt)
      if (rbVal=="Make New Model"){
        makeNewModel(scriptLocation)
      }
      if (rbVal=="Apply Old Model")
        applyOldModel(scriptLocation)
  }
  OK.but <- tkbutton(tt,text="         OK         ",command=OnOK)
  tkgrid(OK.but,columnspan=2)
  tkgrid(tklabel(tt,text="    "))
  tkgrid(tklabel(tt,text="    "))
  tkfocus(tt)
}



