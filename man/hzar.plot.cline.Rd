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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

data

% ## The function is currently defined as
% function (cline, add = FALSE, ...) 
% {
%     if (inherits(cline, "hzar.cline")) 
%         curve(cline$clineFunc(x), add = add, ...)
%     if (inherits(cline, c("hzar.dataGroup", "hzar.fitRequest"))) {
%         dataGroup <- hzar.fit2DataGroup(cline)
%         hzar.plot.obsData(dataGroup, add = add, ...)
%         hzar.plot.cline(hzar.get.ML.cline(dataGroup), add = TRUE, 
%             ...)
%     }
%     if (inherits(cline, c("hzar.obsDataGroup"))) {
%         hzar.plot.obsData(cline, add = add, ...)
%         lapply(cline$data.groups, function(dataGroup) hzar.plot.cline(hzar.get.ML.cline(dataGroup), 
%             add = TRUE, ...))
%     }
%   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
