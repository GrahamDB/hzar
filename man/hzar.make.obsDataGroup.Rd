\name{hzar.make.obsDataGroup}
\alias{hzar.make.obsDataGroup}
\alias{hzar.obsDataGroup}
\title{
  Collect optimizer output from multiple models for analysis.
%%  ~~function to do ... ~~
}
\description{
  This method collects optimizer output based on the same
  \code{\link{hzar.obsData}} object.  It automatically creates
  hzar.dataGroup objects as needed, as well as combining any objects
  which use the same model as determined by
  \code{\link{hzar.sameModel}}.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.obsDataGroup(dataGroups, obsDataGroup = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroups}{
    A list of \code{\link{hzar.dataGroup}} or
    \code{\link{hzar.fitRequest}} objects to include.
    %%     ~~Describe \code{dataGroups} here~~
  }
  \item{obsDataGroup}{
    The \code{hzar.obsDataGroup} to which \code{dataGroups} will be
    added. If \code{NULL}, an empty hzar.obsDataGroup will be added.
    %%     ~~Describe \code{obsDataGroup} here~~
  }
}
% \details{
  
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A hzar.obsDataGroup object.
  \item{data.groups}{A list of hzar.dataGroup objects, each with a unique
    meta models. See \code{\link{hzar.sameModel}}. }
  \item{obsData}{The \code{\link{hzar.obsData}} object shared
    by\code{data.groups}}
  }
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}

%%  ~~who you are~~
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.sameModel}}
  \code{\link{hzar.obsData}}
  \code{\link{hzar.dataGroup}}
}
\examples{
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
mknAdaAmodelFitR$mcmcParam$chainLength <- 5e3;
mknAdaAmodelFitR$mcmcParam$burnin <- 5e2;
mknAdaAmodelFit <- hzar.doFit(mknAdaAmodelFitR)
mknAdaAmodelData <-
  hzar.dataGroup.add(mknAdaAmodelFit);
\dontrun{
mknAdaAmodelFitL <-
    hzar.chain.doSeq(hzar.next.fitRequest(mknAdaAmodelFit), count=3);
}
mknAdaAmodelNull <- hzar.dataGroup.null(mknAdaA);
mknAdaAdGs <- list(clineModel = mknAdaAmodelData,
                   nullModel  = mknAdaAmodelNull);
mknAdaAoDG <- hzar.make.obsDataGroup(mknAdaAmodelFit);
\dontrun{mknAdaAoDG <- hzar.make.obsDataGroup(mknAdaAmodelFitL,mknAdaAoDG);}
mknAdaAoDG <- hzar.make.obsDataGroup(mknAdaAmodelNull,mknAdaAoDG);

mknAdaAoDG <- hzar.copyModelLabels(mknAdaAdGs,mknAdaAoDG);
hzar.plot.cline(mknAdaAoDG);
print(hzar.AICc.hzar.obsDataGroup(mknAdaAoDG));

}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
