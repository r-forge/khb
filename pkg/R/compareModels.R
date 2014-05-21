#' Compute comparable coefficients accross models.
#'
#' The \code{ystandcoef} returns the y-standardized coefficients. This is done by dividing each coefficients by the standard deviation of the latent y dependant variable.
#' 
#' The \code{apecoef} returns the average partial effect of each coefficients. 
#' 
#' @param x A logistic or probit regression. 
#' @return The y-standardized coefficients in which case the \code{"sdy"} attribute of the returned object store the "y" standard deviation or the average partial effects of each covariates.


ystandcoef <- function(x){
	const <- (pi^2)/3
	if(family(x)$link=="probit"){
		const <- 1
	}
	sdy <- sqrt(var(predict(x, type="link"))+const)
	xc <- coef(x)/sdy
	attr(xc, "sdy") <- sdy
	return(xc)
}

#' @rdname ystandcoef

apecoef <- function(x){
	if(family(x)$link=="probit"){
		dform <- dnorm
	}else{
		dform <- dlogis
	}
	return(mean(dform(predict(x, type="link")))*coef(x))
}


getRatioByType <- function(gr, gf, ratio.type, coefFunc=coef){
	cr <- coefFunc(gr)
	cf <-coefFunc(gf)
	uc <- names(cr)[names(cr) %in% names(cf)]
	cr <- cr[uc]
	cf <-cf[uc]
	if(ratio.type=="perc"){
		arc <- 100*(cr-cf)/cr
	} else if (ratio.type=="diff"){
		arc <- cf-cr
	} else{
		arc <- cr/cf
	}
	return(arc)
}

#' @rdname compareModels
#' @param R Integer. If larger than \code{1}, confidence interval of the change are computed using bootstrap. 

coefChange <- function(reduced, full, method="ystand", ratio.type="perc", R=1){
	method <- tolower(method)
	if(method=="ape"){
		coefFunc <- apecoef
	} else if(method=="naive"){
		coefFunc <- coef
	}else if(method=="ystand"){
		coefFunc <- ystandcoef
	} else if(method=="khb"){
		x <- khb(reduced, full)
		reduced <- x$adjusted
		coefFunc <- coef
	}
	
	if(R>1){
		coefFunc <- cmpfun(coefFunc)
		coefChangeInternal <- cmpfun(function(dat, ww){
			ggr <- update(reduced, weights=ww)
			ggf <- update(full, weights=ww)
			return(getRatioByType(ggr, ggf, ratio.type=ratio.type, coefFunc=coefFunc))
		})
		bts <- boot(1:nrow(reduced$model), statistic=coefChangeInternal, sim="ordinary", R=R, stype="f")
		return(bts)
	}

	return(getRatioByType(reduced, full, ratio.type=ratio.type, coefFunc=coefFunc))
}


#' Compare two logistic/probit regression coefficient using different methods.
#' 
#' Compare two logistic/probit regression coefficient using different methods.
#' 
#' The \code{compareModels} function compares the coefficients of two logistic/probit
#' regressions using different methods. \itemize{
#' \item If \code{method="naive"}, the coefficients are compared without any transformation. This method is mainly
#' provided to highlight the effect of other methods. 
#' \item If \code{method="ystand"}, y-standardization is used to compare the coefficients. 
#' \item If \code{method="ape"}, Coefficients change are computed using average partial effects. 
#' \item Finally, if \code{method="khb"}, the KHB method is used to compare the coefficients accross models.
#' }
#' 
#' The \code{coefChange} function only compute the change of coefficients using either
#' the percentage of change or the ratio of change. It can also compute
#' confidence intervals of this change using bootstrapping.
#' 
#' @param reduced The reduced model.
#' @param full The full model
#' @param method Character. The method used to compare coefficients (see
#'   details).
#' @param ratio.type Character specifying how to compare the coefficients. This
#'   should be one of \code{"perc"} (percentage change between the two models),
#'   \code{"ratio"} (coefficient of the reduced model divided the coefficient of the full
#'   one) or \code{"diff"} (difference between the two models).


compareModels <- function(reduced, full, method="ystand", ratio.type="perc"){
	method <- tolower(method)
	if(method=="ape"){
		coefFunc <- apecoef
	} else if(method=="naive"){
		coefFunc <- coef
	}else if(method=="ystand"){
		coefFunc <- ystandcoef
	} else if(method=="khb"){
		x <- khb(reduced, full)
		reduced <- x$adjusted
		coefFunc <- coef
	}
	rtable <- regTable(Reduced=reduced, Full=full, coefFunc=coefFunc, 
				signif.stars=FALSE, se=FALSE, ratio.type=ratio.type, 
				redum="Reduced", fullm="Full")
	print(rtable, quote=FALSE)
}


