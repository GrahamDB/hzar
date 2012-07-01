\name{hzar.plot.obsData}
\alias{hzar.plot.obsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a plot of the observed data points.
}
\description{
  Plots the associated observed frequency versus distance for a variety
  of hzar objects. 
}
\usage{
hzar.plot.obsData(x, type = "p", pch = "+",
                  xlab = "Distance", ylab = hzar.yLabel(x),
                  add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The object from which to extract the observed data to plot.
  }
  \item{type}{
    The plot type for the scatter plot.
  }
  \item{pch}{
    The mark to use to plot the data points.
%%     ~~Describe \code{pch} here~~
  }
  \item{xlab}{
    The x axis label.
    %%     ~~Describe \code{xlab} here~~
  }
  \item{ylab}{
    The y axis label.
    %%     ~~Describe \code{ylab} here~~
  }
  \item{add}{
    Draw on an existing plot.
    %%     ~~Describe \code{add} here~~
  }
  \item{\dots}{
    Additional parameters to pass to plot().
    %%     ~~Describe \code{\dots} here~~
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
%   Software Note
%   %% ~put references to the literature/web site here ~
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
  \code{\link{hzar.extract.obsData}}
}
% \examples{
% ##TODO
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
