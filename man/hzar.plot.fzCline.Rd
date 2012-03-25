\name{hzar.plot.fzCline}
\alias{hzar.plot.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot the 95\% credible cline region for the given locus model.
  %% ~~function to do ... ~~
}
\description{
  Plots the maximum likelihood cline and observed frequency data
  over a the associated fuzzy cline region.  The default region
  is the 95\% credible cline region.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.plot.fzCline(dataGroup, fzCline = hzar.getCredParamRed(dataGroup), type = "p", pch = "+", col = "black", fzCol = "gray", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
    The hzar.dataGroup object for which to generate a fuzzy cline.
    Defaults to a 95\% credible interval region.
%%     ~~Describe \code{dataGroup} here~~
}
\item{fzCline}{
  The hzar.fzCline object to plot.
%%     ~~Describe \code{fzCline} here~~
}
\item{type}{
  The type parameter to pass to hzar.plot.obsData.
%%     ~~Describe \code{type} here~~
}
\item{pch}{
  The plotting character to pass to hzar.plot.obsData.
%%     ~~Describe \code{pch} here~~
}
\item{col}{
  The color to plot the maximum likelihood cline and the observed
  frequencies.
  
%%     ~~Describe \code{col} here~~
}
\item{fzCol}{
  The color to fill the fuzzy cline region with.
%%     ~~Describe \code{fzCol} here~~
}
\item{\dots}{
  Additional parameters to pass to the initial call to plot.
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
  \code{\link{hzar.getCredParamRed}}
  \code{\link{hzar.make.fzCline}}
  \code{\link{plot}}
  \code{\link{hzar.plot.obsData}}
  \code{\link{hzar.plot.cline}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (dataGroup, fzCline = hzar.getCredParamRed(dataGroup), 
%     type = "p", pch = "+", col = "black", fzCol = "gray", ...) 
% {
%     hzar.plot.obsData(dataGroup, col = "transparent", ...)
%     xSeries <- seq(from = par("usr")[1], to = par("usr")[2], 
%         length.out = 109)
%     if (par("xaxs") == "r") 
%         xSeries <- xSeries[2:108]
%     fzCoor <- fzCline$fzCline(xSeries)
%     polygon(x = c(fzCoor$x, rev(fzCoor$x)), y = c(fzCoor$yMin, 
%         rev(fzCoor$yMax)), border = fzCol, col = fzCol)
%     lines(x = xSeries, y = dataGroup$ML.cline$clineFunc(xSeries), 
%         col = col)
%     hzar.plot.obsData(dataGroup, col = col, type = type, pch = pch, 
%         add = TRUE)
%   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
