\name{hzar.getLLCutParam}
\alias{hzar.getLLCutParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get the region of parameter space close to the maximum likelihood
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.getLLCutParam(dataGroups, params, cutValue = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroups}{
%%     ~~Describe \code{dataGroups} here~~
}
  \item{params}{
%%     ~~Describe \code{params} here~~
}
  \item{cutValue}{
%%     ~~Describe \code{cutValue} here~~
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
% function (dataGroups, params, cutValue = 2) 
% {
%     params <- as.character(params)
%     if (inherits(dataGroups, c("hzar.fitRequest", "hzar.dataGroup"))) {
%         dataGroups <- hzar.fit2DataGroup(dataGroups)
%     }
%     else if (is.list(dataGroups)) {
%         if (is.character(names(dataGroups))) {
%             return(do.call(rbind, sapply(names(dataGroups), function(x) hzar.getLLCutParam(dataGroups[[x]], 
%                 params, cutValue), simplify = FALSE)))
%         }
%         return(do.call(rbind, lapply(dataGroups, hzar.getLLCutParam, 
%             params, cutValue)))
%     }
%     else {
%         stop("hzar.getLLCutParam does not understand class of dataGroups.")
%     }
%     data.param = dataGroups$data.param[dataGroups$data.LL$model.LL > 
%         max(dataGroups$data.LL$model.LL - cutValue), ]
%     tempFunc <- function(x) {
%         res <- list(min(data.param[[x]]), max(data.param[[x]]))
%         names(res) <- paste(x, cutValue, "LL", c("Low", "High"), 
%             sep = "")
%         return(res)
%     }
%     return(do.call(data.frame, do.call(c, lapply(params, tempFunc))))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
