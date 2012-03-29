\name{hzar.qScores}
\alias{hzar.qScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate credibility intervals.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.qScores(x, wt, probs = c(0, 0.25, 0.5, 0.75, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{wt}{
%%     ~~Describe \code{wt} here~~
}
  \item{probs}{
%%     ~~Describe \code{probs} here~~
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
% function (x, wt, probs = c(0, 0.25, 0.5, 0.75, 1)) 
% {
%     if (is.data.frame(x)) {
%         res <- list(q = probs)
%         temp.func <- function(index) hzar.qScores(as.numeric(x[, 
%             index]), wt = wt, probs = probs)
%         if (is.character(names(x))) {
%             res2 <- lapply(names(x), temp.func)
%             names(res2) <- names(x)
%             return(do.call(data.frame, c(res, res2)))
%         }
%         return(do.call(data.frame, c(res, lapply(1:(dim(x)[[2]]), 
%             temp.func))))
%     }
%     xS <- sort(unique(x))
%     nS <- length(xS)
%     x0 <- xS[1]
%     xN <- xS[nS]
%     dScale <- (xN - x0)
%     dLen <- xS[2:nS] - xS[1:(nS - 1)]
%     if (is.data.frame(wt)) 
%         wt <- wt[[1]]
%     wt <- wt - max(wt)
%     print("D")
%     scores <- c(0, hzar.qScores.getScores(xSeries = xS, dSeries = dLen/dScale, 
%         raw.x = x, raw.wt = wt))
%     cS <- cumsum(scores/sum(scores))
%     print("E")
%     junk.func <- function(prob, cDSeries, xSeries) {
%         if (length(prob) > 1 || is.list(prob)) 
%             return(sapply(prob, junk.func, cDSeries, xSeries))
%         if (prob < 0) 
%             prob <- 0
%         if (prob > 1) 
%             prob <- 1
%         if (length(which(.Machine$double.eps > (cDSeries - prob)^2)) > 
%             0) 
%             return(mean(xSeries[which(.Machine$double.eps > (cDSeries - 
%                 prob)^2)]))
%         i1 <- rev(which(cDSeries < prob))[[1]]
%         i2 <- which(cDSeries > prob)[[1]]
%         return(xSeries[[i1]] + (prob - cDSeries[[i1]]) * (xSeries[[i2]] - 
%             xSeries[[i1]])/(cDSeries[[i2]] - cDSeries[[i1]]))
%     }
  %   return(sapply(probs, junk.func, cS, xS))
  % }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
