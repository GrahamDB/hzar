\name{hzar.plot.obsData}
\alias{hzar.plot.obsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a plot of the observed data points.
%%  ~~function to do ... ~~
}
\description{
  Plots the associated observed frequency versus distance for a variety
  of hzar objects. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.plot.obsData(x, type = "p", pch = "+", xlab = "Distance", ylab = hzar.yLabel(x), add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The object from which to extract the observed data to plot.
%%     ~~Describe \code{x} here~~
}
\item{type}{
  The plot type for the scatter plot.
%%     ~~Describe \code{type} here~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Graham Derryberry <email: asterion@alum.mit.edu>

}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{plot}}
  \code{\link{ hzar.extract.obsData}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (x, type = "p", pch = "+", xlab = "Distance", ylab = hzar.yLabel(x), 
%     add = FALSE, ...) 
% {
%     oD <- hzar.extract.obsData(x)
%     if (inherits(oD, c("hzar.obsData", "clineSampleData1D", "clineSampleData1DCLT"))) {
%         oDF <- oD$frame
%         yData <- numeric(length(oDF$dist))
%         if (inherits(oD, "clineSampleData1D")) 
%             yData <- oDF$obsFreq
%         if (inherits(oD, "clineSampleData1DCLT")) 
%             yData <- oDF$obsMean
%         if (add) {
%             points(y = yData, x = oDF$dist, type = type, pch = pch, 
%                 ...)
%         }
%         else {
%             plot(y = yData, x = oDF$dist, type = type, pch = pch, 
%                 xlab = xlab, ylab = ylab, ...)
%         }
%     }
%   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
