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

    result1 <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr
    +dmass1+dmass2+ps+pt+py+fancymean+number+median+mean+bscore+sumscore+noa+nda
    +nsa+toa+tda+tsa+percunass+percweakass+percnondirass+onehit+onestronghit
    +onedirecthit+assp+strongassp+dirassp+percp+avgpinten+avgstrongpinten+
    avgdirectpinten+pintenratio+strongpintenratio+dirpintenratio, file1, family=binomial)

### get pvalue matrix from psummary.glm() ###
    object <- psummary.glm(result1)
    test1 <- sort(object[-1,4], decreasing=TRUE)
    
    while (test1[1]>pvalue){
        testnames <- names(test1)
        print(noquote(paste("... removing ",testnames[1])))
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
reducedmodel(filein, 0.05)
#####
#####
#####
filein <- read.csv("pvip_complete.csv")

#Xcorr1 model
#result1 <- glm(v~xcorr1, filein, family=binomial) 

#Xcorrp model
#result1 <- glm(v~xcorrp, filein, family=binomial) 

#Sequest model 
#result1 <- glm(v~xcorr1+charge+mhmass+ionsnum+ionsden+ionsratio, filein, family=binomial) 

#Sequest plus model 
#result1 <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr
+dmass1+dmass2+ps+pt+py, filein, family=binomial)

#Spectral Full Model 
#result1 <- glm(v~xcorrp+charge+dc2+mhmass+ionsnum+ionsden+ionsratio+sp+aa+kr
+dmass1+dmass2+ps+pt+py+fancymean+number+median+mean+bscore+sumscore+noa+nda
+nsa+toa+tda+tsa+percunass+percweakass+percnondirass+onehit+onestronghit
+onedirecthit+assp+strongassp+dirassp+percp+avgpinten+avgstrongpinten+
avgdirectpinten+pintenratio+strongpintenratio+dirpintenratio, filein, family=binomial)

#NEED TO INCORPORATE SPECTRAL REDUCED MODEL
#result1 <- glm(v~xcorrp+noa+percweakass+py+ionsden+strongpintenratio+tda+aa+kr+avgpinten
+avgstrongpinten, filein, family=binomial)

#summary(result1)


test1 <- stepAIC(result1, direction="backward")
print(summary(test1))

test1 <- step(result1, direction="backward")
t2 <- sort(t1[,4], decreasing=T)
