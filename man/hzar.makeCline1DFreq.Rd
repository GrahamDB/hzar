\name{hzar.makeCline1DFreq}
\alias{hzar.makeCline1DFreq}
\alias{hzar.makeCline1DCLT}
\alias{hzar.makeCline1DNormal}
\alias{clineMetaModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a cline model with the requested attributes.
}
\description{
  Constructs a \code{clineMetaModel} object for use with
  \code{\link{hzar.first.fitRequest.old.ML}}.  Said object can be
  further tailored to the specific model desired, or can be used as-is.
}
\usage{
hzar.makeCline1DFreq(data = NULL, scaling = "none", tails = "none",
  direction = NULL) 
hzar.makeCline1DCLT(data = NULL, scaling = "free", tails = "none",
  direction = NULL)
hzar.makeCline1DNormal(data, tails = "none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A \code{\link{hzar.obsData}} object, used to determine cline
    direction and estimate initial values.
    %%     ~~Describe \code{data} here~~
  }
  \item{scaling}{
    Can be one of three strings:
    \itemize{
      \item{"none"}{A model with fixed minimum value 0 and maximum value 1
	is desired.}
      \item{"fixed"}{A model with minimum and maximum values fixed to the
	minimum and maxumimum observed mean values of \code{data} is
	desired.} 
      \item{"free"}{A model with the minimum and maximum value as free
	parameters is desired.}
    }
    %%     ~~Describe \code{scaling} here~~
  }
  \item{tails}{
    Can be one of five strings:
    \itemize{
      \item{"none"}{A model with no exponential tails is desired}
      \item{"right"}{A model with just one exponential tail on the right
	is desired.} 
      \item{"left"}{A model with just one exponential tail on the left is
	desired.} 
      \item{"mirror"}{A model with two exponential tails mirrored about
	the cline center is desired.}
      \item{"both"}{A model with two tails with independent parameters is
	desired.} 
    }
  }
  \item{direction}{
    Can be one of three values:
    \itemize{
      \item{NULL}{Determine direction using \code{data}}
      \item{"ascending"}{A model whose estimates increase as the site
	distance increases is desired.}
      \item{"descending"}{A model whose estimates decrease as the site
	distance increases is desired.}
    }
    %%     ~~Describe \code{direction} here~~
  }
}
\details{
  The \code{clineMetaModel} object returned by
  \code{hzar.makeCline1DNormal} has a slightly diffent structure, due to
  the complexity of the normal cline model.  Use
  \code{\link{hzar.first.fitRequest.gC}} instead of
  \code{\link{hzar.first.fitRequest.old.ML}} to construct the
  \code{\link{hzar.fitRequest}} object needed for \code{\link{hzar.doFit}}.
}
\value{
  A \code{clineMetaModel} object, which is a list with the following
  4 components:
  \item{req}{A boolean function of the model parameters w}
  \item{prior}{Description of 'comp1'}
  \item{func}{Description of 'comp1'}
  \item{parameterTypes}{A list of \code{clineParameter} objects, named
    with the parameter names.
    A \code{clineParameter} object structure:
    \itemize{
      \item{}{Components:}
      \item{val}{The initial or fixed value.}
      \item{w}{The parameter tuning.}
      \item{}{Attributes:}
      \item{"param"}{The parameter name.}
      \item{"fixed"}{TRUE if the parameter is fixed.}
      \item{"limit.lower"}{The parameter minimum finite value.}
      \item{"limit.upper"}{The parameter maximum finite value.}
      \item{"realBTWN01"}{The parameter is restricted to between 0 and
	1.} 
      }
    }
}
\references{
Gay, L., P.-A. Crochet, D. A. Bell, and T. Lenormand. 2008. Comparing clines on molecular and phenotypic traits in hybrid zones: a window on tension zone models. Evolution 62:2789-2806.

Szymura, J., and N. H. Barton. 1986. Genetic analysis of a hybrid zone between the fire-bellied toads, Bombina bombina and B. variegata, near Cracow in souhern Poland. Evolution 40:1141-1159.

Szymura, J., and N. H. Barton. 1991. The genetic structure of the hybrid zone between the fire-bellied toads Bombina bombina and B. variegata: comparisons between transects and between loci. Evolution 45:237-261.

%% ~put references to the literature/web site here ~
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.obsData}} Locality based observed data.
  \code{\link{hzar.doNormalData1DRaw}} Individual observed data.
  \code{\link{hzar.doMorphoSets}} Bernouilli transformed individual observed data.
  
  \code{\link{hzar.model.addBoxReq}} Adding model requirements.
  \code{\link{hzar.first.fitRequest.old.ML}}
  \code{\link{hzar.first.fitRequest.gC}}
}
\examples{
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
mknAdaAm <-
  hzar.makeCline1DFreq(mknAdaA, scaling="fixed",tails="none");
str(mknAdaAm);
## compile using this method:
mknAdaAFR <-
hzar.first.fitRequest.old.ML(model=mknAdaAm, obsData=mknAdaA)
hzar.plot.cline(mknAdaAFR)

data(manakinMorphological)
data(manakinLocations)
mknCCbt <-
  hzar.doMorphoSets("collar.color",
                    tDist=manakinLocations,
                    tDLocCol="LocalityID",
                    tDDistCol="distance",
                    tValues=manakinMorphological,
                    tVLocCol="Locality")$collar.color
mknCCbtM <-
  hzar.makeCline1DFreq(mknCCbt, scaling="fixed",tails="none");
str(mknCCbt);
## compile using this method:
mknCCbtFR <-
  hzar.first.fitRequest.old.ML(model=mknCCbtM, obsData=mknCCbt)
hzar.plot.cline(mknCCbtFR)

data(manakinMorphological);
data(manakinLocations);
mknBLg <-
  hzar.doNormalData1DRaw(
    hzar.mapSiteDist(siteID=manakinLocations$LocalityID,
                     distance=manakinLocations$distance),
    traitSite=manakinMorphological$Locality,
    traitValue=manakinMorphological$beard.length)
hzar.plot.obsData(mknBLg)
mknBLgM <-
  hzar.makeCline1DNormal(mknBLg, tails="none");

## Emulate scaling="fixed"

hzar.meta.fix(mknBLgM)$muL <- TRUE
hzar.meta.fix(mknBLgM)$muR <- TRUE
## Use variance of site at the left (site A) and right (site L). 
hzar.meta.init(mknBLgM)$varL<-mknBLg$frame["A","var"]
hzar.meta.init(mknBLgM)$varR<-mknBLg$frame["L","var"]
hzar.meta.fix(mknBLgM)$varL <- TRUE
hzar.meta.fix(mknBLgM)$varR <- TRUE

## compile using this method:
mknBLgFR <- hzar.first.fitRequest.gC(gModel=mknBLgM, obsData=mknBLg)
hzar.plot.cline(mknBLgFR)

data(manakinMorphological);
data(manakinLocations);
bcD=manakinLocations[manakinMorphological$Locality,
                     "distance"]
bcV=manakinMorphological$belly.color
bcV=bcV[!is.na(bcD)];bcD=bcD[!is.na(bcD)];
bcD=bcD[!is.na(bcV)];bcV=bcV[!is.na(bcV)];
mknBCclt <-
  hzar.doCLTData1DRaw( distance=bcD,traitValue=bcV)
\dontrun{
hzar.plot.obsData(mknBCclt)
}
mknBCcltM <-
  hzar.makeCline1DCLT(mknBCclt, scaling="free", tails="none");


}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
