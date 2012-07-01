\name{hzar.cov.rect}
\alias{hzar.cov.rect}
\alias{hzar.cov.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a covariance matrix for the cline optimizer.
}
\description{
  These methods are intended to generate covariance matrices suitable
  for use with MCMCmetrop1R.
}
\usage{
hzar.cov.rect(clineLLfunc, param.lower, param.upper,
              pDiv = 11, random = 0, passCenter = FALSE)
hzar.cov.mcmc(clineLLfunc, mcmcRaw,
              pDiv = 15, random = 10000, passCenter = FALSE)
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
  }  
  \item{mcmcRaw}{
    A mcmc object used to refine the covariance matrix.
  }
  
}
\details{
  This method is adaptive, refining the pDiv and random parameters until
  either it can generate a useable matrix without too high of a memory
  cost, or random > 1e9 (use a 1 billion or more samples).
}
\value{
  A square matrix with a width equal to the number of free parameters.  
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
% \examples{
% ##
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
