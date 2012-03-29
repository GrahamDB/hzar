\name{hzar.chain.doSeq}
\alias{hzar.chain.doSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Repeatedly run the optimizer in series to tune the optimizer.
%%  ~~function to do ... ~~
}
\description{
  To be written.
}
\usage{
hzar.chain.doSeq(hzar.request, count = 3, collapse = FALSE, announce.complete = "Chain Complete")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hzar.request}{
%%     ~~Describe \code{hzar.request} here~~
}
  \item{count}{
%%     ~~Describe \code{count} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
}
  \item{announce.complete}{
%%     ~~Describe \code{announce.complete} here~~
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
Graham Derryberry \email{asterion@alum.mit.edu}
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

% ## The function is currently defined as
% function (hzar.request, count = 3, collapse = FALSE, announce.complete = "Chain Complete") 
% {
%     if (collapse) {
%         mcmcParam = hzar.request$mcmcParam
%         mcmcParam$chainLength = mcmcParam$chainLength * count
%         mcmcParam$burnin = mcmcParam$burnin * count
%         mdlParam = hzar.request$modelParam
%     }
%     hzar.results <- list()
%     for (iter in 1:count) {
%         print(iter)
%         if (iter > 1) {
%             hzar.request <- NULL
%             print(summary(try(hzar.request <- hzar.next.fitRequest(hzar.results[[iter - 
%                 1]]))))
%         }
%         if (!inherits(hzar.request, c("hzar.fitRequest"))) {
%             warning("Failed to generate next fit request. Returning successful runs.")
%             return(hzar.results)
%         }
%         hzar.results[[iter]] <- hzar.request
%         print(summary(try(hzar.results[[iter]] <- hzar.doFit(hzar.request))))
%     }
%     print(announce.complete)
%     if (collapse) {
%         raw.data <- do.call(rbind, lapply(hzar.results, function(x) x$mcmcRaw))
%         rawMCMC <- mcmc(data = raw.data, thin = thin(hzar.results[[count]]$mcmcRaw), 
%             start = mcmcParam$burnin + 1)
%         names(rawMCMC) <- names(hzar.results[[count]]$mcmcRaw)
%         return(list(hzar.make.fitRequest(mdlParam, hzar.results[[count]]$cM, 
%             hzar.results[[count]]$llFunc, mcmcParam, mcmcRaw = rawMCMC, 
%             TRUE, prod(as.logical(lapply(test.fits, attr, "fit.success"))))))
%     }
%     return(hzar.results)
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
