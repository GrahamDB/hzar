\name{hzar.cov.rect}
\alias{hzar.cov.rect}
\alias{hzar.cov.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a covariance matrix for the cline optimizer.
%%  ~~function to do ... ~~
}
\description{
  These methods are intended to generate covariance matrices suitable
  for use with MCMCmetrop1R.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.cov.rect(clineLLfunc, param.lower, param.upper, pDiv = 11, random = 0, passCenter = FALSE)
hzar.cov.mcmc(clineLLfunc, mcmcRaw, pDiv = 15, random = 10000, passCenter = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clineLLfunc}{
    The log likelihood function of the parameters.
    %%     ~~Describe \code{clineLLfunc} here~~
  }
  \item{param.lower}{
    The minimum boundary of the region of parameter space to consider.
    %%     ~~Describe \code{param.lower} here~~
  }
  \item{param.upper}{
    The maximum boundary of the region of parameter space to consider.
    %%     ~~Describe \code{param.upper} here~~
  }
  \item{pDiv}{
    If generating a covariance matrix using a lattice, the lattice
    should have this many points on an edge. 
    %%     ~~Describe \code{pDiv} here~~
  }
  \item{random}{
    Use \code{random} number of points drawn from a uniform likelihood
    space to generate the covariance matrix.  If 0, use a lattice to
    generate the covariance matrix.
    %%     ~~Describe \code{random} here~~
  }
  \item{passCenter}{
    Should weighted mean of the parameter space be returned.
    %%     ~~Describe \code{passCenter} here~~
  }  
  \item{mcmcRaw}{
    A mcmc object used to refine the covariance matrix.
    %%     ~~Describe \code{mcmcRaw} here~~
  }
  
}
\details{
  This method is adaptive, refining the pDiv and random parameters until
  either it can generate a useable matrix without too high of a memory
  cost, or random > 1e9 (use a 1 billion or more samples).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A square matrix with a width equal to the number of free parameters.  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{MCMCmetrop1R}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (clineLLfunc, param.lower, param.upper, pDiv = 11, random = 0, 
%     passCenter = FALSE) 
% {
%     if (random > 1e+09) {
%         stop("Covariance matrix calculation with random sampling requested with far too many samples.  Stopping.")
%     }
%     if (random > 0) {
%         data.mat <- list(dTheta = prod(abs(as.numeric(param.upper) - 
%             as.numeric(param.lower)))/random, data = hzar.gen.rParam.uniform(param.lower, 
%             param.upper, random))
%     }
%     else {
%         if ((pDiv^length(param.lower)) > 1e+06) {
%             warning("Covariance matrix calculation requested for too complex of a lattice structure. Switching to random sampling using ten thousand samples.")
%             return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
%                 random = 10000, passCenter = passCenter))
%         }
%         data.mat <- hzar.gen.samples.rect(param.lower, param.upper, 
%             pDiv)
%     }
%     param.names <- names(data.mat$data)
%     data.wt <- hzar.eval.clineLL(data.mat$data, clineLLfunc)
%     data.mat$data <- data.mat$data[data.wt > -1e+06, ]
%     data.wt <- data.wt[data.wt > -1e+06]
%     MIN.DATA <- (1 + length(param.upper))
%     if (length(data.wt) < MIN.DATA) {
%         if (random > 0) {
%             return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
%                 random = 2 * random, passCenter = passCenter))
%         }
%         else {
%             return(hzar.cov.rect(clineLLfunc, param.lower, param.upper, 
%                 pDiv = pDiv + 1, passCenter = passCenter))
%         }
%     }
%     while (sum(data.wt > -723) < MIN.DATA) {
%         if (sum(data.wt > 609) > 0) {
%             if (random > 0) {
%                 return(hzar.cov.rect(clineLLfunc, param.lower, 
%                   param.upper, random = 2 * random, passCenter = passCenter))
%             }
%             else {
%                 return(hzar.cov.rect(clineLLfunc, param.lower, 
%                   param.upper, pDiv = pDiv + 1, passCenter = passCenter))
%             }
%         }
%         data.wt <- data.wt + 100
%     }
%     VDATA <- cov.wt(x = cbind(data.mat$data, model.LL = data.wt), 
%         wt = exp(data.wt))
%     VMATRIX <- VDATA$cov
%     mat.scaled <- VMATRIX
%     if (passCenter) 
%         return(list(cov = mat.scaled[param.names, param.names], 
%             center = VDATA$center[param.names]))
%     return(mat.scaled[param.names, param.names])
%   }
% ## The function is currently defined as
% function (clineLLfunc, mcmcRaw, pDiv = 15, random = 10000, passCenter = FALSE) 
% {
%     mcmc.nm <- colnames(mcmcRaw)
%     pL <- lapply(mcmc.nm, function(x) {
%         min(mcmcRaw[, x])
%     })
%     names(pL) <- mcmc.nm
%     pU <- lapply(mcmc.nm, function(x) {
%         max(mcmcRaw[, x])
%     })
%     names(pU) <- mcmc.nm
%     return(hzar.cov.rect(clineLLfunc, pL, pU, pDiv, random, passCenter))
%   }

}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
