\name{hzar.yLabel}
\alias{hzar.yLabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Suggest a y axis label
  %%  ~~function to do ... ~~
}
\description{
  Suggests a y axis label based on the hzar object passed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.yLabel(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
% function (x) 
% {
%     oD <- hzar.extract.obsData(x)
%     if (inherits(oD, "clineSampleData1D")) {
%         return("Frequency")
%     }
%     if (inherits(oD, "clineSampleData1DCLT")) {
%         return("Trait Value")
%     }
%     if (is.list(oD)) {
%         if (prod(as.logical(lapply(oD, inherits, what = "clineSampleData1D"))) == 
%             1) 
%             return("Frequency")
%         if (prod(as.logical(lapply(oD, inherits, what = "clineSampleData1DCLT"))) == 
%             1) 
%             return("Trait Value")
%     }
%     return("")
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
