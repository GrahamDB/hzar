\name{hzar.profile.dataGroup}
\alias{hzar.profile.dataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a likelihood profile for a single parameter
}
\description{
  TODO
}
\usage{
hzar.profile.dataGroup(dG, parameter, pVals = NULL, pDivs = NULL, nDiv =
20, appeture = NULL, doPar = FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dG}{
    A \code{\link{hzar.dataGroup}} of the fitted model
  }
  \item{parameter}{
    The name of the parameter to profile.
    %%     ~~Describe \code{parameter} here~~
  }
  \item{pVals}{
    The parameter values at which to calculate the profile.
    Automatically sorted and reduced to unique values.  If NULL,
    default to the midpoints of ranges specified by \code{pDiv}.
    %%     ~~Describe \code{pVals} here~~
  }
  \item{pDivs}{
    A list of ranges of parameter values from the supplied
    \code{\link{hzar.dataGroup}} to use for initializing the profile
    fit. If numeric, e.g. \code{c(0,1,2,4)}, treat as if a list of pairs
    of consequetive values, e.g. \code{list(c(0,1),c(1,2),c(2,4))}. If
    NULL, generate using \code{pVals} and \code{appeture}.  If
    \code{pVals} is also NULL, generate using \code{\link{equal.count}},
    with \code{n = nDiv}.

    %%     ~~Describe \code{pDivs} here~~
  }
  \item{nDiv}{
    The number of points at which to calculate the profile. Ignored if
    either \code{pDivs} or \code{pVals} are not NULL.
    %%     ~~Describe \code{nDiv} here~~
  }
  \item{appeture}{
    The maximum parameter delta from \code{pVals} to use for
    initializing the profile fit from the supplied
    \code{\link{hzar.dataGroup}}. Ignored if \code{pDivs} is not NULL or
    \code{pVals} is not \code{numeric}. If used and NULL, default to the
    median deviation between consecutive values of \code{pVals}. 
    %%     ~~Describe \code{appeture} here~~
  }
  \item{doPar}{
    Should \code{\link{\%dopar\%}} be used to construct fit requests in
    parallel? 
    %%     ~~Describe \code{doPar} here~~
  }
  \item{\dots}{
    Arguments to pass to \code{\link{hzar.multiFitRequest}}.
  }
}
\details{
  This method does not actually do the fitting of the likelihood
  profile. See the example for a useable workflow.

  
}
\value{
  A list of \code{\link{hzar.fitRequest}} objects to be fitted.
}
\references{
  There is probably a reference to how to do this somewhere...
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
\seealso{
  \code{\link{hzar.multiFitRequest}}
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.doFit.multi}}
  \code{\link{hzar.make.obsDataGroup}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
%% % Add one or more standard keywords, see file 'KEYWORDS' in the
%% % R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
