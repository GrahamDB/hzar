\name{hzar-package}
\alias{hzar-package}
\alias{hzar}
\docType{package}
\title{
  Hybrid Zone Analysis using R
}
\description{
  A collection of tools for modeling the shape of 1 dimensional clines.
}
\details{
\tabular{ll}{
Package: \tab hzar\cr
Type: \tab Package\cr
Version: \tab 0.1-7\cr
Date: \tab 2012-03-29\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
 Graham Derryberry
Maintainer:  Graham Derryberry <asterion@alum.mit.edu>
%~~ The author and/or maintainer of the package ~~
}
\references{
  Brumfield, R. T., R. W. Jernigan, D. B. McDonald, and
  M. J. Braun. 2001. Evolutionary implications of divergent clines in an
  avian (Manacus: Aves) hybrid zone. Evolution 55:2070-2087. 
  
  Gay, L., P.-A. Crochet, D. A. Bell, and T. Lenormand. 2008. Comparing
  clines on molecular and phenotypic traits in hybrid zones: a window on
  tension zone models. Evolution 62:2789-2806.
  
  Szymura, J., and N. H. Barton. 1986. Genetic analysis of a hybrid zone
  between the fire-bellied toads, Bombina bombina and B. variegata, near
  Cracow in souhern Poland. Evolution 40:1141-1159.  
  
  Szymura, J., and N. H. Barton. 1991. The genetic structure of the
  hybrid zone between the fire-bellied toads Bombina bombina and
  B. variegata: comparisons between transects and between
  loci. Evolution 45:237-261. 
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
% \seealso{
%   \code{\link[MCMCpack:MCMCpack-package]{MCMCpack}}
% }
\examples{
%~~ simple examples of the most important functions ~~
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
hzar.plot.obsData(mknAdaA);
mknAdaAmodel <-
  hzar.makeCline1DFreq(mknAdaA, scaling="fixed",tails="none");
mknAdaAmodel <-
  hzar.model.addBoxReq(mknAdaAmodel,-30,600);
mknAdaAmodelFitR <-
   hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,
                                mknAdaA,
                                verbose=FALSE);
mknAdaAmodelFitR$mcmcParam$chainLength <- 2e3;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e2;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
plot(hzar.mcmc.bindLL(mknAdaAmodelFit))
mknAdaAmodelData <-
  hzar.dataGroup.add(mknAdaAmodelFit);
\dontrun{
mknAdaAmodelData <-
  hzar.dataGroup.add(
    mknAdaAmodelData,
    hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit)));
hzar.plot.cline(mknAdaAmodelData);
hzar.plot.fzCline(mknAdaAmodelData);
}
print(hzar.getLLCutParam(mknAdaAmodelData,c("center","width")));
mknAdaAmodelNull <- hzar.dataGroup.null(mknAdaA);
mknAdaAdGs <- list(clineModel = mknAdaAmodelData,
                   nullModel  = mknAdaAmodelNull);
mknAdaAoDG <- hzar.make.obsDataGroup(mknAdaAdGs);
mknAdaAoDG <- hzar.copyModelLabels(mknAdaAdGs,mknAdaAoDG);
hzar.plot.cline(mknAdaAoDG);
print(hzar.AICc.hzar.obsDataGroup(mknAdaAoDG));

}
