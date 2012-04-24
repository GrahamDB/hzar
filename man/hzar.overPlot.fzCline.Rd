\name{hzar.overPlot.fzCline}
\alias{hzar.overPlot.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot a cline region for multiple models and / or loci. 
}
\description{
  To be written.
}
\usage{
hzar.overPlot.fzCline(dataGroupSet,
                      fzClineSet = sapply(dataGroupSet,
                                          hzar.getCredParamRed,
                                          simplify = FALSE),
                      type = "p",
                      fzDens = 8,
                      fzShadeAngle =
                        ((1:length(dataGroupSet)) * 180)
                         \%/\% (1 + length(dataGroupSet)),
                      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroupSet}{
    The list of \code{\link{hzar.dataGroup}} objects to generate
    \code{\link{hzar.fzCline}} objects for plotting.
  }
  \item{fzClineSet}{
    The list of \code{\link{hzar.fzCline}} objects to plot.
  }
  \item{type}{
    Passed to \code{\link{hzar.plot.obsData}}.
  }
  \item{fzDens}{
    Density of the shading lines. See \code{\link{polygon}}.
  }
  \item{fzShadeAngle}{
    Angle of the shading lines. See \code{\link{polygon}}.
  }
  \item{\dots}{
    Additional parameters to pass to \code{\link{hzar.plot.fzCline}}. 
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
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
  \code{\link{hzar.getCredParamRed}}
  \code{\link{hzar.plot.fzCline}}
}
\examples{
##TOdO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
