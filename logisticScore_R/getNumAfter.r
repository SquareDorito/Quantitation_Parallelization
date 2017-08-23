getNumAfter<-function(string,afterThis){ ## may need mods for cterm & nterm
if (afterThis=="STY*"){

  BreakOnParen=strsplit(string,"(",fixed=T)
  index=pmatch(afterThis,BreakOnParen[[1]]) #returns NA if no match -- hence if/else stmt
  BreakOnParen=strsplit(BreakOnParen[[1]][index],")",fixed=T)

  index=pmatch(afterThis,BreakOnParen[[1]]) 
  PsudoNumber=strsplit(BreakOnParen[[1]][index],afterThis,fixed=T)

  len=length(PsudoNumber[[1]])
  if(len>1){
    number=PsudoNumber[[1]][len]
  }
  else{
    number=0
  }

}
else if(afterThis=='Cterm-pep='||afterThis=='Nterm-pep='){
  breakOn=strsplit(string,afterThis,fixed=T)
  breakOnSpace=strsplit(breakOn[[1]][length(breakOn[[1]])],' ',fixed=T)
  PsudoNumber=breakOnSpace[[1]][1]
  if(is.na(as.numeric(PsudoNumber))){number=0}
  else{number=PsudoNumber}

}
else{
  BreakOnSpaces=strsplit(string," ",fixed=T)
  index=pmatch(afterThis,BreakOnSpaces[[1]]) #returns NA if no match -- hence if/else stmt 
  PsudoNumber=strsplit(BreakOnSpaces[[1]][index],afterThis,fixed=T)
  len=length(PsudoNumber[[1]])
  if(len>1){
    number=PsudoNumber[[1]][len]
  }
  else{
    number=0
  }
}
return(as.numeric(number))

}