\name{hzar.doCLTData1DRaw}
\alias{hzar.doCLTData1DRaw}
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
    The distance of the population site of the sample.  Samples at the
    same distance are treated as being from the same population site.
  }
  \item{traitValue}{
    The value of the trait of the indidual sampled.
  }
}
\details{
  If for any population site, there is only a small number of samples
  taken, warnings will be issued.

  If at any population site, the sample variance is 0, a warning is
  issued, and additional variance is included by estimating the amount
  of variance ignored due to measurement error.

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
