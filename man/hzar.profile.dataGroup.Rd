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
    %%     ~~Describe \code{parameter} here~~
  }
  \item{pVals}{
    %%     ~~Describe \code{pVals} here~~
  }
  \item{pDivs}{
    %%     ~~Describe \code{pDivs} here~~
  }
  \item{nDiv}{
    %%     ~~Describe \code{nDiv} here~~
  }
  \item{appeture}{
    %%     ~~Describe \code{appeture} here~~
  }
  \item{doPar}{
    %%     ~~Describe \code{doPar} here~~
  }
  \item{\dots}{
    Arguments to pass to \code{\link{hzar.multiFitRequest}}
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
