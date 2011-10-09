\name{fitClineModel}
\alias{fitClineModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitClineModel(model, sampleData, verbose = 10000, burnin = 1000, mcmc = 1e+06, thin = 100, seed = list(NA, 1), rejectLL = -1e+08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{sampleData}{
%%     ~~Describe \code{sampleData} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{mcmc}{
%%     ~~Describe \code{mcmc} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{rejectLL}{
%%     ~~Describe \code{rejectLL} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, sampleData, verbose = 10000, burnin = 1000, 
    mcmc = 1e+06, thin = 100, seed = list(NA, 1), rejectLL = -1e+08) 
{
    myClineFrame = list(data = sampleData)
    attr(myClineFrame, "mcmc.count") <- 0
    if (identical(is.null(seed), TRUE)) {
        seed = list(NA, attr(clineFrame, "mcmc.count") + 1)
    }
    if (is.list(seed) && length(seed) > 1) 
        attr(myClineFrame, "mcmc.count") <- seed[[2]]
    myParams <- splitParameters(model$parameterTypes)
    oldFormals <- formals(model$func)
    tttFormals <- oldFormals[names(myParams$init)]
    names(tttFormals) <- names(myParams$init)
    formals(model$func) <- c(tttFormals, myParams$fixed)
    formals(model$req) <- c(tttFormals, myParams$fixed)
    formals(model$prior) <- c(tttFormals, myParams$fixed)
    myRejectionLL <- rejectLL
    clineLLFunc <- function(theta, meta.model, obsData) {
        if (!do.call(meta.model$req, as.list(theta))) 
            return(myRejectionLL)
        model = do.call(meta.model$func, as.list(theta))
        thetaLL = do.call(meta.model$prior, as.list(theta))
        result <- obsData$model.LL(model) + thetaLL
        if (identical(is.finite(result), TRUE)) 
            return(result)
        return(myRejectionLL)
    }
    VMATRIX <- NULL
    pMinS.cssp <- NULL
    pMaxS.cssp <- NULL
    if (sum(names(myParams$fixed) \%in\% "pMin") > 0) 
        pMinS.cssp <- myParams$fixed$pMin
    if (sum(names(myParams$fixed) \%in\% "pMax") > 0) 
        pMaxS.cssp <- myParams$fixed$pMax
    sampleModels <- cssp(sampleData$frame, attr(model, "tails"), 
        pMinS = NULL, pMaxS = NULL)[, names(myParams$init)]
    if (dim(sampleModels)[[1]] > 2) {
        sampleModels <- subset(sampleModels, do.call(model$req, 
            as.list(sampleModels[names(myParams$init)])))
        for (iter in seq(dim(sampleModels)[[1]])) {
            sampleModels[iter, "model.LL"] <- clineLLFunc(sampleModels[iter, 
                names(myParams$init)], model, sampleData)
        }
        try(VMATRIX <- cov(sampleModels))
        if (!identical(is.null(VMATRIX), TRUE)) {
            print(VMATRIX)
            bounded.param <- lapply(model$parameterTypes[names(myParams$init)], 
                attr, "realBTWN01") == TRUE
            counter.inv <- diag(ifelse(c(bounded.param, FALSE), 
                sign(VMATRIX["model.LL", ]), 1/(VMATRIX["model.LL", 
                  ])))
            dimnames(counter.inv) <- list(rownames(VMATRIX), 
                colnames(VMATRIX))
            VMATRIX <- (counter.inv \%*\% VMATRIX \%*\% counter.inv)
            VMATRIX <- VMATRIX[names(myParams$init), names(myParams$init)]
        }
    }
    print("C")
    print(VMATRIX)
    myClineFrame$mcmc <- MCMCmetrop1R(fun = clineLLFunc, logfun = "TRUE", 
        burnin = burnin, mcmc = mcmc, thin = thin, theta.init = myParams$init, 
        tune = myParams$tune, meta.model = model, force.samp = TRUE, 
        obsData = sampleData, seed = seed, V = VMATRIX, verbose = verbose, 
        optim.control = list(fnscale = -1, trace = 0, REPORT = 10, 
            maxit = 5000), optim.method = "L-BFGS-B", optim.lower = myParams$lower, 
        optim.upper = myParams$upper)
    oldParams <- model$parameterTypes
    colnames(myClineFrame$mcmc) <- names(myParams$init)
    model.info = as.data.frame(myClineFrame$mcmc)
    nSamples = length(model.info[[names(myParams$init)[[1]]]])
    if (length(myParams$fixed > 0)) {
        for (iter in names(myParams$fixed)) {
            model.info[[iter]] = rep(myParams$fixed[[iter]], 
                nSamples)
        }
    }
    for (iter in 1:nSamples) {
        model.info[iter, "model.LL"] <- clineLLFunc(model.info[iter, 
            names(myParams$init)], model, sampleData)
    }
    myClineFrame$model = model
    myClineFrame$allClines = model.info
    myClineFrame$maxLL = max(model.info$model.LL)
    myClineFrame$maxLL.theta = subset(model.info, model.info$model.LL == 
        max(model.info$model.LL))[1, ]
    myClineFrame$maxLL.theta <- myClineFrame$maxLL.theta[names(formals(model$func))]
    myClineFrame$maxLL.cline = do.call(model$func, as.list(myClineFrame$maxLL.theta))
    class(myClineFrame) <- "clineModelFrame"
    return(myClineFrame)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
