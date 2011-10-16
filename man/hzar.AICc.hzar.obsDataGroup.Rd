\name{hzar.AICc.hzar.obsDataGroup}
\alias{hzar.AICc.hzar.obsDataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the corrected AIC score (AICc) table for the
  hzar.obsDataGroup collection.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.AICc.hzar.obsDataGroup(obsDataGroup, label = "AICc", show.count = FALSE, show.param = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsDataGroup}{
%%     ~~Describe \code{obsDataGroup} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{show.count}{
%%     ~~Describe \code{show.count} here~~
}
  \item{show.param}{
%%     ~~Describe \code{show.param} here~~
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

## The function is currently defined as
function (obsDataGroup, label = "AICc", show.count = FALSE, show.param = FALSE) 
{
    model.clines <- lapply(obsDataGroup$data.groups, hzar.get.ML.cline)
    scores <- as.numeric(lapply(model.clines, hzar.AICc.hzar.cline, 
        sum(obsDataGroup$obsData$frame$n)))
    result <- data.frame(AIC = scores)
    if (!identical(is.null(label), TRUE)) {
        names(result) <- label
        colnames(result) <- label
    }
    if (!identical(is.null(names(obsDataGroup$data.groups)), 
        TRUE)) {
        rownames(result) <- names(obsDataGroup$data.groups)
    }
    if (show.count) {
        count = as.numeric(lapply(model.clines, function(x) length(x$param.free)))
        result <- cbind(result, count = count)
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
