\name{hzar.make.fzCline}
\alias{hzar.make.fzCline}
\alias{hzar.fzCline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a \code{hzar.fzCline} object
}
\description{
  Compile a list of \code{\link{hzar.cline}} objects into a
  \code{hzar.fzCline} object. 
}
\usage{
hzar.make.fzCline(clineList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clineList}{
    A list of clines.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{hzar.fzCline} object.
  \item{clines}{\code{clineList}}
  \item{listFuncInt}{ A function of a scalar x}
  \item{fzCline}{ A function of a numeric series over which
    \code{clineList} is evaluated and the series of maximum and minimum
    values is returned. }
}
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

% \seealso{
%   \code{\link{hzar.fzCline}} 
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% % ## The function is currently defined as
% % function (clineList) 
% % {
% %     res <- list(clines = clineList, listFuncInt = function(xVal, 
% %         funcList = res$clines) {
% %         yList <- as.numeric(lapply(funcList, function(x, u) x$clineFunc(u), 
% %             u = xVal))
% %         return(data.frame(x = xVal, yMin = min(yList), yMax = max(yList)))
% %     }, fzCline = function(xVals, listFunc = res$listFuncInt) {
% %         xVals <- as.numeric(xVals)
% %         if (length(xVals > 1)) return(do.call(rbind, lapply(xVals, 
% %             listFunc)))
% %         if (length(xVals < 1)) return(numeric(0))
% %         return(listFunc(xVals))
% %     })
% %     return(res)
% %   }
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
