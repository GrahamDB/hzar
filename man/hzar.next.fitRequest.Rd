\name{hzar.next.fitRequest}
\alias{hzar.next.fitRequest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a new fitRequest using data from another fitRequest.
%%  ~~function to do ... ~~
}
\description{
  The method is the glue for parallel runs of sequential chains.
  It returns a ready to run \code{\link{hzar.fitRequest}} object based
  on the \code{\link{hzar.fitRequest}} supplied.  If
  \code{oldFitRequest} had already been successfully run, this method's
  result will be dependent on the prior run.  If not, this method's
  result will be a request for an independent run.
}
\usage{
hzar.next.fitRequest(oldFitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oldFitRequest}{
    A \code{\link{hzar.fitRequest}} object.
  }
}
\details{
  If \code{oldFitRequest} describes a successful run, a new covariance
  matrix is generated, modelParam$init is updated to the covariance
  matrix center, and the lecuyer seed channel incremented by 1. If the
  mersenne twister was used previously, the lecuyer random number
  generator is requested on channel 2.

  If \code{oldFitRequest} does not describe a successful run, everything
  is copied, except the lecuyer seed channel incremented by 10 to
  prevent overlapping. If the mersenne twister was used previously, the
  lecuyer random number generator is requested on channel 11. 
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{\link{hzar.fitRequest}} object.
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
  Graham Derryberry \email{asterion@alum.mit.edu}

}
% \note{
%   The comment about overlapping seed channels is a result of the
%   following scenario:

%   \code{ 
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.chain.doSeq}}
  \code{\link{hzar.cov.mcmc}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (oldFitRequest) 
% {
%     seedChannel <- 1
%     if (is.list(oldFitRequest$mcmcParam$seed)) {
%         seedChannel = oldFitRequest$mcmcParam$seed[[2]]
%     }
%     if (identical(attr(oldFitRequest, "fit.run"), TRUE)) {
%         seedChannel <- seedChannel + 1
%     }
%     else {
%         seedChannel <- seedChannel + 10
%     }
%     mcmcParam = hzar.make.mcmcParam(oldFitRequest$mcmcParam$chainLength, 
%         oldFitRequest$mcmcParam$burnin, oldFitRequest$mcmcParam$verbosity, 
%         oldFitRequest$mcmcParam$thin, seedChannel)
%     mdlParam <- oldFitRequest$modelParam
%     covMatrix <- oldFitRequest$cM
%     if (identical(attr(oldFitRequest, "fit.success"), TRUE)) {
%         mcmcSubset <- oldFitRequest$mcmcRaw[sample(dim(oldFitRequest$mcmcRaw)[[1]]), 
%             ]
%         subLL <- hzar.eval.clineLL(mcmcSubset, oldFitRequest$llFunc)
%         covData <- hzar.cov.mcmc(oldFitRequest$llFunc, mcmcSubset[subLL > 
%             max(subLL - 4), ], passCenter = TRUE)
%         covMatrix <- covData$cov
%         new.center <- covData$center[names(mdlParam$init)]
%         if (oldFitRequest$llFunc(new.center) > 1e-06) 
%             mdlParam$init <- new.center
%     }
%     return(hzar.make.fitRequest(mdlParam, covMatrix, oldFitRequest$llFunc, 
%         mcmcParam))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
