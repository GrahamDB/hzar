\name{hzar.getCredCut}
\alias{hzar.getCredCut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get the minimum log likelihood of sufficient credibility.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.getCredCut(dataGroup, rejectionPercent = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
  \item{rejectionPercent}{
%%     ~~Describe \code{rejectionPercent} here~~
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
% function (dataGroup, rejectionPercent = 0.05) 
% {
%     if (inherits(dataGroup, c("hzar.fitRequest", "hzar.dataGroup"))) {
%         dataGroup <- hzar.fit2DataGroup(dataGroup)
%     }
%     else if (is.list(dataGroup) && length(dataGroup) > 0 && inherits(dataGroup[[1]], 
%         c("hzar.fitRequest", "hzar.dataGroup"))) {
%         warning("Only calculating LL cut for first element in list")
%         dataGroup <- hzar.fit2DataGroup(dataGroups[[1]])
%     }
%     model.relLL = exp(sort(dataGroup$data.LL$model.LL - dataGroup$ML.cline$logLike))
%     credibleLLspace <- data.frame(LL = sort(dataGroup$data.LL$model.LL), 
%         percentile = cumsum(model.relLL/sum(model.relLL)))
%     credible.LLcut <- min(subset(credibleLLspace, credibleLLspace$percentile > 
%         rejectionPercent)$LL)
%     return(credible.LLcut)
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
