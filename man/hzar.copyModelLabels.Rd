\name{hzar.copyModelLabels}
\alias{hzar.copyModelLabels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Copy names from one hzar object to another.
}
\description{
  Set the \code{\link{names}} of the list of
  \code{\link{hzar.dataGroup}} objects contained in a
  \code{\link{hzar.obsDataGroup}} object using the names from either a
  named list of \code{\link{hzar.dataGroup}} objects or another
  \code{\link{hzar.obsDataGroup}} object.
}
\usage{
hzar.copyModelLabels(group1, group2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{group1}{
    An object from which to extract model labels.
  }
  \item{group2}{
    A \code{\link{hzar.obsDataGroup}} into which the model labels need
    to be inserted.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  The updated \code{group2}.
}

\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.obsDataGroup}}
  \code{\link{names}}
}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% % ## The function is currently defined as
% % function (group1, group2) 
% % {
% %     if (inherits(group1, c("hzar.obsDataGroup"))) 
% %         group1 <- group1$data.groups
% %     if (inherits(group2, c("hzar.obsDataGroup"))) {
% %         group2$data.groups <- hzar.copyModelLabels(group1, group2$data.groups)
% %         return(group2)
% %     }
% %     if (!is.list(group1) || !is.list(group2)) 
% %         stop("group1 or group2 of malformed type")
% %     names(group2) <- as.character(lapply(group2, function(x) names(group1)[which(as.logical(lapply(group1, 
% %         hzar.sameModel, x)))[1]]))
% %     return(group2)
% %   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
