\name{hzar.plot.fzCline}
\alias{hzar.plot.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the 95\% credible cline region for the given locus model.
}
\description{
  Plots the maximum likelihood cline and observed frequency data
  over a the associated fuzzy cline region.  The default region
  is the 95\% credible cline region.
}
\usage{
hzar.plot.fzCline(dataGroup,
                  fzCline = hzar.getCredParamRed(dataGroup),
                  type = "p", pch = "+",
                  col = "black", fzCol = "gray", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
    The hzar.dataGroup object for which to generate a fuzzy cline.
    Defaults to a 95\% credible interval region.
  }
  \item{fzCline}{
    The hzar.fzCline object to plot.
  }
  \item{type}{
    The type parameter to pass to hzar.plot.obsData.
  }
  \item{pch}{
    The plotting character to pass to hzar.plot.obsData.
  }
  \item{col}{
    The color to plot the maximum likelihood cline and the observed
    frequencies.
  }
  \item{fzCol}{
    The color to fill the fuzzy cline region with.
  }
  \item{\dots}{
    Additional parameters to pass to the initial call to plot.
    
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
  \code{\link{hzar.make.fzCline}}
  \code{\link{plot}}
  \code{\link{hzar.plot.obsData}}
  \code{\link{hzar.plot.cline}}
}
\examples{
##TODO
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
