\name{hzar.gen.cline}
\alias{hzar.gen.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a hzar.cline object using the given parameters and model.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.gen.cline(free.parameters, fitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{free.parameters}{
%%     ~~Describe \code{free.parameters} here~~
}
  \item{fitRequest}{
%%     ~~Describe \code{fitRequest} here~~
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
% function (free.parameters, fitRequest) 
% {
%     context <- hzar.extract.modelPrep(fitRequest)
%     cline.func <- NULL
%     cline.param <- c(as.list(free.parameters), context$param.fixed)
%     if (do.call(context$model.req, as.list(free.parameters))) {
%         cline.func <- do.call(context$model.gen, as.list(free.parameters))
%     }
%     return(hzar.make.cline(free.parameters, parameters = cline.param, 
%         func = cline.func, LL = context$llFunc(free.parameters)))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
