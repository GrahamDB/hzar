\name{hzar.model.addBoxReq}
\alias{hzar.model.addBoxReq}
\alias{hzar.model.addCenterRange}
\alias{hzar.model.addMaxCenter}
\alias{hzar.model.addMinCenter}
\alias{hzar.model.addMaxDelta}
\alias{hzar.model.addMaxWidth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add parameter restriction clauses to cline model
  %%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.model.addBoxReq(meta.model, low, high)
hzar.model.addCenterRange(meta.model, low, high)
hzar.model.addMaxCenter(meta.model, maxValue)
hzar.model.addMinCenter(meta.model, minValue)
hzar.model.addMaxDelta(meta.model, maxValue)
hzar.model.addMaxWidth(meta.model, maxValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meta.model}{
%%     ~~Describe \code{meta.model} here~~
}
  \item{low}{
%%     ~~Describe \code{low} here~~
}
  \item{high}{
%%     ~~Describe \code{high} here~~
}
  \item{minValue}{
%%     ~~Describe \code{high} here~~
}
  \item{maxValue}{
%%     ~~Describe \code{high} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Graham Derryberry \email{asterion@alum.mit.edu}

%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (meta.model, low, high) 
% {
%     return(hzar.model.addCenterRange(hzar.model.addMaxDelta(hzar.model.addMaxWidth(meta.model, 
%         high - low), high - low), low, high))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
