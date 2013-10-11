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
  
  If for any locality, there is only a small number of samples
  taken, warnings will be issued.

  If at any locality, the sample variance is 0, a warning is
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
  \code{\link{hzar.obsData}} For a description of this object structure.
  \code{\link{hzar.makeCline1DNormal}} Make models of normal data.
  \code{\link{hzar.first.fitRequest.gC}} For compiling those models.
  \code{\link{hzar.doFit}} For fitting the compiled models.
  \code{\link{hzar.doFit.multi}} For quickly fitting multiple models.
  \code{\link{hzar.mcmc.bindLL}} For viewing the mcmc trace.
}
\examples{
data(manakinMorphological);
data(manakinLocations);
mknBL <-
  hzar.doNormalData1DRaw(
    hzar.mapSiteDist(siteID=manakinLocations$LocalityID,
                     distance=manakinLocations$distance),
    traitSite=manakinMorphological$Locality,
    traitValue=manakinMorphological$beard.length)
hzar.plot.obsData(mknBL)
mknBLm <-
  hzar.makeCline1DNormal(mknBL, tails="none");
mknBLm <-
  hzar.model.addBoxReq(mknBLm,-30,600);
## A quick way to reduce the number of variables

## Assume the observed means to left and right are the
## population means
hzar.meta.fix(mknBLm)$muL <- TRUE
hzar.meta.fix(mknBLm)$muR <- TRUE

## Make the initial variance to the left and right match the
## observed local variance at the left (site A) and right (site L). 
hzar.meta.init(mknBLm)$varL<-mknBL$frame["A","var"]
hzar.meta.init(mknBLm)$varR<-mknBL$frame["L","var"]

## Assume the observed variance to left and right is the
## the respective population variance. 
hzar.meta.fix(mknBLm)$varL <- TRUE
hzar.meta.fix(mknBLm)$varR <- TRUE

## Now mknBLm has only 3 free parameters instead of 7.

mknBLFR <-
   hzar.first.fitRequest.gC(gModel=mknBLm,
                            obsData=mknBL,
                            verbose=FALSE);
mknBLFR$mcmcParam$chainLength <- 2e3;
mknBLFR$mcmcParam$burnin <- 5e2;
mknBLF <- hzar.doFit(mknBLFR)
plot(hzar.mcmc.bindLL(mknBLF))
\dontrun{
mknBLFR2 <- hzar.next.fitRequest(mknBLF)
## Do more fitting
mknBLFR2$mcmcParam$chainLength <- 1e5;
mknBLFR2$mcmcParam$burnin <- 1e3;
mknBLF2 <- hzar.doFit(mknBLFR2)
plot(hzar.mcmc.bindLL(mknBLF2))
}
}

% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ observed }  
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
