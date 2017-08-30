
#### FUNCTIONS 
psummary.glm <- function (object) 
{
    # revision of summary.glm to give matrix of Pvalues to 
    # recreate stepwise regression in Stata, To avoid reduction using AIC values
    est.disp <- FALSE
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        aliased <- is.na(coef(object))
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coef.p, Inf)
            dimnames(coef.table) <- list(names(coef.p), dn)
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        aliased <- is.na(coef(object))
        df.f <- length(aliased)
    }
    ans <- data.matrix(coef.table)    
    return(ans)
}

reducedmodel <- function (file1, pvalue) 
{
### rewriting a step function based on P(Z) pvalue instead of AIC values (to replicate stata)
### initial logistic regression full spectral model
### returns glm.object of stepwise logistic reduction model below the given pvalue

    result1 <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr+dmass1+dmass2+ps+pt+py+fancymean+number+median+mean+bscore+sumscore+noa+nda+nsa+toa+tda+tsa+percunass+percweakass+percnondirass+onehit+onestronghit+onedirecthit, file1, family=binomial)

### get pvalue matrix from psummary.glm() ###
    object <- psummary.glm(result1)
    test1 <- sort(object[-1,4], decreasing=TRUE)
    
    while (test1[1]>pvalue){
        testnames <- names(test1)
        #print(testnames)
        print(noquote(paste("... removing ",testnames[1])))
	    #print("---------------------------------")
        #print(summary(result1))
        funcout<-"v~"
        for (i in 2:length(testnames)){
            if (i == 2)
               funcout = paste(funcout,testnames[i],sep="")
            else
		    funcout = paste(funcout,"+",testnames[i],sep="")
    	  }
	      result1 <- glm(as.formula(funcout), file1, family=binomial)
	      object <- psummary.glm(result1)
        test1 <- sort(object[-1,4], decreasing=TRUE)	
    }
   
    return(result1)
}

predict <- function (coeff, file1) {
### function for creating the prediction based on the model generated coefficients
### returns vector of predictions

	p1 = vector(mode="numeric",length=nrow(file1)); 
	p2 = p1;
	for (i in 2:length(coeff)) {
		namecol <- names(coeff[i])
		if (!is.na(coeff[i])){
			p1 = p1 + (coeff[i]*file1[,namecol])
		}
	}
	p1 <- p1+coeff[1]
	#print(coeff[1])
	#print(p1[1])
	p2 <- exp(p1)/(1+exp(p1))
	p2[p1>=100]=1
	#print(p2[1])
	
	return(p2)
}
logreg2<-function(curdir,filein){
#### MAIN SORT OF ####
### FILE INPUT FOR NOW ####
#curdir =  "//proteome.biomed.brown.edu/User Files/Anthony/Anthony/final_combinations/fnished matrix/"

library("ROC")

## sets R's current working directory
setwd(curdir)
#filein <- read.csv("bsa_complete2.csv")
#filein <- read.csv("mcp5_complete2.csv")
#filein <- read.csv("pvip_complete.csv")

#Sequest model 
sequest <- glm(v~xcorr1+charge+mhmass+ionsnum+ionsden+ionsratio+sp, filein, family=binomial)
#print(curdir)
xc <- glm(v~xcorr1, filein, family=binomial)
xcp <- glm(v~xcorrp, filein, family=binomial)
save(sequest,file=paste(curdir,"/sequest.RData",sep=""))
#Sequest plus model 
sequestplus <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr
+dmass1+dmass2+ps+pt+py, filein, family=binomial)
save(sequestplus,file=paste(curdir,"/sequestplus.RData",sep=""))
#Spectral Full Model 
full <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr+dmass1+dmass2+ps+pt+py+fancymean+number+median+mean+bscore+sumscore+noa+nda+nsa+toa+tda+tsa+percunass+percweakass+percnondirass+onehit+onestronghit+onedirecthit, filein, family=binomial)
save(full,file=paste(curdir,"/full.RData",sep=""))
#SPECTRAL REDUCED MODEL
reduced <- reducedmodel(filein, 0.05)
save(reduced,file=paste(curdir,"/reduced.RData",sep=""))
## Summary of coefficients and P(z)
#summary(result1)

### GETS PREDICTION
sequestout <- predict(coefficients(sequest),filein)
sequestplusout <- predict(coefficients(sequestplus),filein)
fullout <- predict(coefficients(full),filein)
reducedout <- predict(coefficients(reduced),filein)
xcorrout <- predict(coefficients(xc),filein)
xcorrpout <- predict(coefficients(xcp),filein)

### NEED TO LOAD ROC LIBRARY
### TODO::::USE ROCR library, more diagnostic plots!
### produces AUC and ROC plots

rxcorr <- rocdemo.sca(filein[,1], filein[,"xcorr1"], dxrule.sca)
rxcorrp <- rocdemo.sca(filein[,1], filein[,"xcorrp"], dxrule.sca)
rsequest <- rocdemo.sca(filein[,1], sequestout, dxrule.sca)
rsequestplus <- rocdemo.sca(filein[,1], sequestplusout, dxrule.sca)
rfull <- rocdemo.sca(filein[,1], fullout, dxrule.sca)
rreduce <- rocdemo.sca(filein[,1], reducedout, dxrule.sca)

fprXcorr=getFDR(filein[,"rev_database"],xcorrout)
fprXcorrp=getFDR(filein[,"rev_database"],xcorrpout)
fprSequest=getFDR(filein[,"rev_database"],sequestout)
fprSequestP=getFDR(filein[,"rev_database"],sequestplusout)
fprSpec=getFDR(filein[,"rev_database"],fullout)
fprRed=getFDR(filein[,"rev_database"],reducedout)

par(mfrow=c(3,4))
#print(rxcorr)
#print(fprXcorr)
#print(filein[,1])
#print(xcorrout)
#print(coefficients(xc))

plot(rxcorr,line=TRUE)
title(main=paste("xcorr:",round(AUC(rxcorr),4)))

plot(rxcorrp,line=TRUE)
title(main=paste("xcorrp:",round(AUC(rxcorrp),4)))

plot(rsequest,line=TRUE)
title(main=paste("sequest:",round(AUC(rsequest),4)))

plot(rsequestplus,line=TRUE)
title(main=paste("sequestplus:",round(AUC(rsequestplus),4)))

plot(rfull,line=TRUE)
title(main=paste("spectral full:",round(AUC(rfull),4)))

plot(rreduce,line=TRUE)
title(main=paste("spectral reduced:",round(AUC(rreduce),4)))

plot(fprXcorr[,1],fprXcorr[,2],"l")
title(main="FPR: Xcorr")

plot(fprXcorrp[,1],fprXcorrp[,2],'l')
title(main="FPR: Xcorr'")

plot(fprSequest[,1],fprSequest[,2],'l')
title(main="FPR: Sequest")

plot(fprSequestP[,1],fprSequestP[,2],'l')
title(main="FPR: Sequest Plus")

plot(fprSpec[,1],fprSpec[,2],'l')
title(main="FPR: Spectral Full")

plot(fprRed[,1],fprRed[,2],'l')
title(main="FPR: Spectral Reduced")

#gives area underneath the roc curve
#AUC(r1)
}

