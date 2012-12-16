\name{hzar.doCLTData1DRaw}
\alias{hzar.doCLTData1DRaw}
\alias{hzar.doNormalData1DRaw}
\alias{hzar.mapSiteDist}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a \code{hzar.obsData} object using a table of individual traits.
}
\description{
  Create a \code{\link{hzar.obsData}} object using a table of individual
  traits. 
}
\usage{
hzar.doCLTData1DRaw(distance, traitValue)
hzar.doNormalData1DRaw(site.dist, traitSite, traitValue)
hzar.mapSiteDist(siteID, distance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
    The distance of the sampling site.  For
    \code{hzar.doCLTData1DRaw}, samples at the same distance are treated
    as being from the same sampling site.
  }
  \item{traitValue}{
    The value of the trait of the individual sampled.
  }
  \item{traitSite}{
    The id of site where the individual was found.
  }
  \item{site.dist}{
    A named vector mapping site id codes to the distance of the sampling
    site. The function \code{hzar.mapSiteDist} returns a suitable vector.
  }
  \item{siteID}{
    The list of id codes associated with the sampling site.  This list
    should be identical in length to \code{distance}, each entry must be
    unique, and the order of the sites referenced must be identical for
    \code{distance} and \code{siteID}.
  }
}
\details{
  For \code{hzar.doCLTData1DRaw}:
  
  If for any population site, there is only a small number of samples
  taken, warnings will be issued.

  If at any population site, the sample variance is 0, a warning is
  issued, and additional variance is included by estimating the amount
  of variance ignored due to measurement error.

  For \code{hzar.doNormalData1DRaw}:

  Use the helper function \code{hzar.mapSiteDist} to generate
  \code{site.dist}.

  The \code{\link{hzar.obsData}} object created is meant for use with
  the models constructed using \code{\link{hzar.makeCline1DNormal}}.
}
\value{
  A \code{\link{hzar.obsData}} object, using the site dinstances and
  sample means and variances as calculated from the values given.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.obsData}}
}
% \examples{
% ##TODO
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ observed }  
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
