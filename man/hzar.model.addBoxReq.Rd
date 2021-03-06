\name{hzar.model.addBoxReq}
\alias{hzar.model.addBoxReq}
\alias{hzar.model.addCenterRange}
\alias{hzar.model.addMaxCenter}
\alias{hzar.model.addMinCenter}
\alias{hzar.model.addMaxDelta}
\alias{hzar.model.addMaxWidth}
\alias{hzar.model.addMaxVariance}
\alias{hzar.model.addMuRange}
\alias{hzar.model.addNormalBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add parameter restriction clauses to cline model
  %%  ~~function to do ... ~~
}
\description{
  The intent of these methods is to assist the optimizer in exploring
  the model parameter space by instructing it to ignore models that are
  not interesting.  For example, if all of the sampled localities are
  in a region 100km wide, then a cline width of 110km is probably not
  interesting.  A cline width of 500km in that scenario would definitely
  not be interesting at all.
}
\usage{
hzar.model.addBoxReq(meta.model, low, high)
hzar.model.addCenterRange(meta.model, low, high)
hzar.model.addMaxCenter(meta.model, maxValue)
hzar.model.addMinCenter(meta.model, minValue)
hzar.model.addMaxDelta(meta.model, maxValue)
hzar.model.addMaxWidth(meta.model, maxValue)
hzar.model.addMaxVariance(meta.model, maxValue)
hzar.model.addNormalBox(meta.model, left, right, bottom,top)
hzar.model.addMuRange(meta.model, low, high)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta.model}{
    The \code{\link{clineMetaModel}} object to modify.
  }
  \item{minValue}{
    The smallest value to consider.
  }
  \item{maxValue}{
    The greatest value to consider.
  }
  \item{left}{
    The leftmost location to consider.
  }
  \item{right}{
    The rightmost location to consider.
  }
  \item{bottom}{
    The least trait value to consider.
  }
  \item{top}{
    The greatest trait value to consider.
  }
  The following arguments specifier a range in distances, with ascending
  values from left to right. 
  \item{low}{
    The leftmost location to consider.
  }
  \item{high}{
    The rightmost location to consider.
  }
}
\details{
  The three center methods only add requirements to the center
  parameter. Likewise, \code{hzar.model.addMaxWidth} only adds
  a maximum width requirement. In constrast, \code{hzar.model.addMaxDelta}
  adds a maximum value requirement to any and all delta parameters
  present in \code{meta.model}.

  \code{hzar.model.addBoxReq} adds requirements to any and all of
  the parameters center, width, deltaM, deltaL, and deltaR.  The center
  requirements are the same as calling
  \code{hzar.model.addCenterRange(meta.model, low, high)}. The remaining
  parameters are required to have a maximum value of \code{high-low}.
  
}
\value{
  The modified \code{\link{clineMetaModel}} object.
}

\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{clineMetaModel}}
}
\examples{
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
mknAdaAmodel <-
  hzar.makeCline1DFreq(mknAdaA, scaling="fixed",tails="none");
mknAdaAmodelB <-
  hzar.model.addBoxReq(mknAdaAmodel,-30,600);
mknAdaAmodel$req;
mknAdaAmodelB$req;
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
