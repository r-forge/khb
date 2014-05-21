
regTable <- function(..., digits=max(3L, getOption("digits") - 2L), signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, se=TRUE, coefFunc=coef, ratio.type=NULL, redum=1, fullm=2, khb.resid=TRUE){
	sigfunc <- function(pv){
		return(symnum(pv, corr = FALSE, na = FALSE, 
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                  symbols = c("***", "**", "*", ".", " ")))
	}
	roundedDigit <- function(val){
		val <- zapsmall(val, digits)
		acs <- abs(val)
		if (any(ia <- is.finite(acs))) {
			digmin <- 1 + if (length(acs <- acs[ia & acs != 0])) 
				floor(log10(range(acs[acs != 0], finite = TRUE)))
			else 0
			xx <- format(round(val, max(1L, digits - 
				digmin)), digits = digits)
		}else{
			xx <- format(val, digits = digits)
		}
		return(xx)
	}
	coeffSignif <- function(rl){
		if(signif.stars){
			allcof <- summary(rl)$coefficients
			cof <- allcof[, 1]
			significance <- allcof[, 4]
			mycoeff <- paste(
					roundedDigit(cof),
					sigfunc(significance)
					)
			names(mycoeff) <- names(cof)
		}
		else{
			mycoeff <- roundedDigit(coefFunc(rl))
		}
		if(se){
			nc <- length(mycoeff)
			withse <- character(2*nc)
			withse[(1:nc)*2-1] <- mycoeff
			names(withse)[(1:nc)*2-1] <- names(mycoeff)
			withse[(1:nc)*2] <- paste0("(", roundedDigit(allcof[,2]), ")")
			names(withse)[(1:nc)*2] <- paste("SE_", names(mycoeff), sep="")
			return(withse)
		}else{
			return(mycoeff)
		}
		
		
	}
	allreg <- list(...)
	allregraw <- allreg
	if(is.null(names(allreg))){
		names(allreg) <- paste("(", 1:length(allreg), ")", sep="")
	}
	allregcoef <- list()
	diffcoef <- NULL
	for(rr in names(allreg)){
		allregcoef[[rr]] <- coeffSignif(allreg[[rr]])
		diffcoef <- c(diffcoef, names(allregcoef[[rr]]))
	}
	diffcoef <- unique(diffcoef)
	regmat <- matrix("", nrow=length(diffcoef), ncol=length(allreg)+!is.null(ratio.type))
	rownames(regmat) <- diffcoef
	##colnames(regmat) <- names(allreg)
	for(rr in 1:length(allreg)){
		arc <- allregcoef[[names(allreg)[rr]]]
		regmat[names(arc), rr] <- arc
	}
	if(!is.null(ratio.type)){
		arc <- getRatioByType(allreg[[redum]], allreg[[fullm]], ratio.type=ratio.type, coefFunc=coefFunc)
		regmat[names(arc), length(allreg)+1] <- roundedDigit(arc)
	}
	if(se){
		rownames(regmat)[(1:(nrow(regmat)/2))*2] <- ""
	}
	
	rowcond <- grepl("KHB__KHB__", rownames(regmat), fixed=TRUE)
	rownames(regmat)[rowcond] <- gsub("\\`?KHB__KHB__(.*)\\`?$", "Resid(\\1)", rownames(regmat)[rowcond])
	if(!khb.resid){
		regmat <- regmat[!rowcond, ]
	}
	qualmat <- matrix("", ncol=length(allreg), nrow=7)
	rownames(qualmat) <- c("Pseudo R2", "Dev.", "Null", "Chisq", "Sig", "Dl", "BIC")
	qualmat["Pseudo R2", ] <- sapply(allreg, function(x) roundedDigit(NagelkerkeR2(x)$R2))
	qualmat["Dev.", ] <- sapply(allreg, function(x) roundedDigit(deviance(x)))
	qualmat["Null", ] <- sapply(allreg, function(x) roundedDigit(x$null.deviance))
	qualmat["Chisq", ] <- sapply(allreg, function(x) {sx <- summary(x); roundedDigit(sx$null.deviance - sx$deviance)})
	qualmat["Sig", ] <- sapply(allreg, function(x) {sx <- summary(x);
											dl <- sx$df.null - sx$df.residual;
											sigfunc(pchisq(sx$null.deviance - sx$deviance, dl, lower.tail=FALSE))
											}
	)
	qualmat["Dl", ] <- sapply(allreg, function(x) {sx <- summary(x);  paste(sx$df.null - sx$df.residual, sep="")})
	qualmat["BIC", ] <- sapply(allreg, function(x) {sx <- summary(x); roundedDigit(sx$deviance+(sx$df.null - sx$df.residual)*log(nobs(x)))})
	
	finalmat <- rbind(rep("", length(allreg)), qualmat)
	if(!is.null(ratio.type)) finalmat <- cbind(finalmat, rep("", nrow(finalmat)))
	finalmat <- rbind(regmat, finalmat)
	colnames(finalmat) <- c(names(allreg), ratio.type)
	return(finalmat)
}

