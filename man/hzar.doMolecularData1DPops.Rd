\name{hzar.doMolecularData1DPops}
\alias{hzar.doMolecularData1DPops}
\alias{hzar.obsData}
\alias{hzar.doCLTData1DPops}
\alias{hzar.doNormalData1DPops}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a hzar.obsData object using summary data about each locality
}
\description{
  Creates a hzar.obsData object using the observations given.  The
  likelihood function used is chosen based on the method called.
}
\usage{
hzar.doMolecularData1DPops(distance, pObs, nEff,
                           siteID=paste("P",1:length(distance),sep=""),
                           ylim=extendrange(c(0,1)))
hzar.doCLTData1DPops(distance, muObs, varObs, nEff)
hzar.doNormalData1DPops(distance, muObs, varObs, nEff,
                        siteID=paste("P",1:length(distance),sep=""),
                        ylim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  All arguments should be of the same length.
  \item{distance}{
    The distance of each locality. If the same distance is given
    multiple times, then multiple localities are assumed to be at that
    distance.
  }
  \item{pObs}{
    The observed frequency at each locality.
  }
  \item{nEff}{
    The effective number of samples observed at each locality.
  }
  \item{muObs}{
    The mean trait value observed at each site.
  }
  \item{varObs}{
    The trait variance observed at each site.
  }
  \item{ylim}{
    The ylim to use when plotting the observed data.
  }
  \item{siteID}{
    The identifier to use for each sampling site.
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
  \item{frame}{A data.frame composed of the arguments.}
  \item{model.LL}{A function of one argument that returns a log
    likelihood.  The argument is a function of distance that estimates
    either frequency or trait value as appropriate.}
}
\references{
  Szymura, J., and N. H. Barton. 1986. Genetic analysis of a hybrid zone
  between the fire-bellied toads, Bombina bombina and B. variegata, near
  Cracow in souhern Poland. Evolution 40:1141-1159. 

  Szymura, J., and N. H. Barton. 1991. The genetic structure of the
  hybrid zone between the fire-bellied toads Bombina bombina and
  B. variegata: comparisons between transects and between
  loci. Evolution 45:237-261. 
}
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{manakinMolecular}}
  \code{\link{hzar.plot.obsData}}
}
\examples{
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
print(str(mknAdaA));
hzar.plot.obsData(mknAdaA);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
