
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils::packageDescription("khb"))
	if(utils::packageVersion("khb")$minor %% 2 == 0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	if(!is.null(descr$Built)){
		builtDate <- paste(" (Built: ", strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2], ")", sep="")
	}else{
		builtDate <- ""
	}
	packageStartupMessage("This is khb ", state, " version ", descr$Version, builtDate)
	packageStartupMessage("\nIf you use the KHB library, please cite it using:")
	##x <- citation("khb")
	sapply(strwrap("Matthias Studer (2014). khb: Comparing nonlinear regression models. R package version 0.1.", exdent=4), packageStartupMessage)
	packageStartupMessage("\nAnd if you use the KHB method:")
	x <-   "Karlson KB, Holm A and Breen R (2012). Comparing Regression Coefficients Between Same-sample Nested Models Using Logit and Probit: A New Method. Sociological Methodology, 42(1), pp 286-313."
	sapply(strwrap(x, exdent=4), packageStartupMessage)
	invisible()
}
