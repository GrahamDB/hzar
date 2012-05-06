\name{hzar.getCredParam}
\alias{hzar.getCredParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Select a credible subset of a collection of parameters.
}
\description{
  Select the subset of the generated samples wuth a likelihood greater
  that the result of \code{\link{hzar.getCredCut}}. 
}
\usage{
hzar.getCredParam(dataGroup, rejectionPercent = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
    The \code{\link{hzar.dataGroup}} object to analyze.
  }
  \item{rejectionPercent}{
    The proportion of the cumalitive likelihood distribuition to
    reject. 
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{data.frame}} of the subset of the generated samples. 
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

\seealso{
  \code{\link{data.frame}}
  \code{\link{hzar.getCredCut}}
  \code{\link{hzar.dataGroup}}
}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% % ## The function is currently defined as
% % function (dataGroup, rejectionPercent = 0.05) 
% % {
% %     if (inherits(dataGroup, c("hzar.fitRequest", "hzar.dataGroup"))) {
% %         dataGroup <- hzar.fit2DataGroup(dataGroup)
% %     }
% %     else if (is.list(dataGroup) && length(dataGroup) > 0 && inherits(dataGroup[[1]], 
% %         c("hzar.fitRequest", "hzar.dataGroup"))) {
% %         warning("Only calculating LL cut for first element in list")
% %         dataGroup <- hzar.fit2DataGroup(dataGroup[[1]])
% %     }
% %     credible.LLcut <- hzar.getCredCut(dataGroup, rejectionPercent)
% %     return(dataGroup$data.param[dataGroup$data.LL$model.LL > 
% %         credible.LLcut, ])
% %   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
