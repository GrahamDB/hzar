\name{hzar.doMolecularData1DPops}
\alias{hzar.doMolecularData1DPops}
\alias{hzar.obsData}
\alias{hzar.doCLTData1DPops}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a hzar.obsData object using summary data about each population.
%%  ~~function to do ... ~~
}
\description{
  Creates a hzar.obsData object using the observations given.  The
  likelihood function used is chosen based on the method called.
  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.doMolecularData1DPops(distance, pObs, nSamples)
hzar.doCLTData1DPops(distance, muObs, varObs, nSamples)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  All arguments should be of the same length.
  \item{distance}{
    The distance of each population. If the same distance is given
    multiple times, then multiple sites are assumed to be at that
    distance.
  }
  \item{pObs}{
    The observed frequency at each site.
  }
  \item{nSamples}{
    The numbered of samples observed at each site.
  }
  \item{muObs}{
    The mean trait value observed at each site.
  }
  \item{varObs}{
    The trait variance observed at each site.
  }
}
\details{
  For \code{hzar.doCLTData1DPops}, \code{varObs} must not be less than
  zero, and should be greater than zero. If equal to zero, the method
  will attempt to estimate the number of significant digits in the
  observed trait value, and use that to calculate additional variance
  due to measurement error at each site, and add that variance to the
  observed variance for each site.
}
\value{
  A \code{hzar.obsData} object with the following structure.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  \item{frame}{A data.frame composed of the arguments.}
  \item{model.LL}{A function of one argument that returns a log
    likelihood.  The argument is a function of distance that estimates
    either frequency or trait value as appropriate.}
  %% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
% function (distance, pObs, nSamples) 
% {
%     if ((length(distance) != length(pObs)) || (length(distance) != 
%         length(nSamples)) || (length(pObs) != length(nSamples))) {
%         stop("Distance, pObs and nSamples are not of the same length.")
%     }
%     obj <- list(frame = data.frame(dist = distance, obsFreq = pObs, 
%         n = nSamples))
%     obj$model.LL <- function(model.func) {
%         pEst = model.func(obj$frame$dist)
%         res <- sampleLikelihoodMolecularPop(pEst = as.numeric(pEst), 
%             pObs = as.numeric(obj$frame$obsFreq), N = as.numeric(obj$frame$n))
%         result <- sum(res)
%         if (is.na(result)) 
%             return(-1e+08)
%         return(result)
%     }
%     class(obj) <- c("clineSampleData1D", "hzar.obsData")
%     return(obj)
%   }
% function (distance, muObs, varObs, nSamples) 
% {
%     if ((length(distance) != length(muObs)) || (length(distance) != 
%         length(nSamples)) || (length(distance) != length(varObs))) {
%         stop("Distance, muObs, varObs and nSamples are not all of the same length.")
%     }
%     if (sum(nSamples < 2) > 0) 
%         stop("There must be at least two samples per population!")
%     if (sum(nSamples < 10) > 0) 
%         warning("Some populations have less than 10 samples.  More samples would be nice.")
%     if (sum(varObs == 0) > 0) {
%         warning("Some of the population samples have a variance of 0.  Adding in estimated variance due to measurement error.")
%         m.err <- 5/3 * 10^-quantile(1 + as.numeric(lapply(muObs, 
%             function(nV) {
%                 min((-1:12)[round(nV, digits = -1:12) == nV])
%             })), probs = 0.75)[[1]]
%         varObs <- varObs + m.err * m.err
%     }
%     obj <- list(frame = data.frame(dist = distance, obsMean = muObs, 
%         obsVariance = varObs, n = nSamples))
%     obj$model.LL <- function(model.func) {
%         muEst = model.func(obj$frame$dist)
%         res <- logLikeCLTEstPopMean(muEst = as.numeric(muEst), 
%             muObs = as.numeric(obj$frame$obsMean), varObs = as.numeric(obj$frame$obsVariance), 
%             nObs = as.numeric(obj$frame$n))
%         result <- sum(res)
%         if (is.na(result)) 
%             return(-1e+08)
%         return(result)
%     }
%     class(obj) <- c("clineSampleData1DCLT", "hzar.obsData")
%     return(obj)
%   }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
