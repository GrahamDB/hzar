\name{hzar.plot.cline}
\alias{hzar.plot.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a plot of the cline.
%%  ~~function to do ... ~~
}
\description{
  Plots a line representing the expected frequency versus distance for
  the given object.  For hzar.dataGroup and hzar.obsDataGroup objects,
  plots the observed data backing the model. For hzar.obsDataGroup
  objects, plots the maximum likelihood cline for each model.
}
\usage{
hzar.plot.cline(cline, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cline}{
    A hzar.cline, hzar.dataGroup or hzar.obsDataGroup object.
}
  \item{add}{
    Add to an existing plot if TRUE.
}
\item{\dots}{
  Arguments to be passed to methods, such as graphical
  parameters (see \code{\link{plot}}).
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}


% \references{
% %% ~put references to the literature/web site here ~
% }

\author{
  Graham Derryberry \email{asterion@alum.mit.edu}

}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{plot}}
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
mknAdaAcline <- hzar.gen.cline(list(center=300,width=10),
                               mknAdaAmodelFitR);
hzar.plot.cline(mknAdaAmodelFitR);
hzar.plot.cline(mknAdaAcline,add=TRUE);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
