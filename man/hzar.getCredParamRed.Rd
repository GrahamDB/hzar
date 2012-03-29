\name{hzar.getCredParamRed}
\alias{hzar.getCredParamRed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a 95\% credibility hzar.fzCline object
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.getCredParamRed(dataGroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
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
  \code{\link{hzar.make.fzCline}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (dataGroup) 
% {
%     junk <- hzar.getCredParam(dataGroup)
%     junk1 <- PP.fzCline.reduce(junk)
%     junk1 <- PP.fzCline.reduce(junk1)
%     res <- hzar.make.fzCline(lapply(1:(dim(junk1)[[1]]), function(x) hzar.gen.cline(junk1[x, 
%         ], dataGroup)))
%     return(res)
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
