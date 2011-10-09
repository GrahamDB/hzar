\name{hzar.make.clineLLfunc.old.bayes}
\alias{hzar.make.clineLLfunc.old.bayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.clineLLfunc.old.bayes(param.free.names, param.fixed, param.check.func, meta.cline.func, model.LL, prior.LL, LLrejectedModel = -1e+08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param.free.names}{
%%     ~~Describe \code{param.free.names} here~~
}
  \item{param.fixed}{
%%     ~~Describe \code{param.fixed} here~~
}
  \item{param.check.func}{
%%     ~~Describe \code{param.check.func} here~~
}
  \item{meta.cline.func}{
%%     ~~Describe \code{meta.cline.func} here~~
}
  \item{model.LL}{
%%     ~~Describe \code{model.LL} here~~
}
  \item{prior.LL}{
%%     ~~Describe \code{prior.LL} here~~
}
  \item{LLrejectedModel}{
%%     ~~Describe \code{LLrejectedModel} here~~
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
function (param.free.names, param.fixed, param.check.func, meta.cline.func, 
    model.LL, prior.LL, LLrejectedModel = -1e+08) 
{
    model.req <- param.check.func
    model.gen <- meta.cline.func
    eval.priorLL <- prior.LL
    eval.clineLL <- model.LL
    myRejectionLL <- LLrejectedModel
    old.formals <- formals(model.gen)
    if (length(old.formals) != (length(param.free.names) + length(param.fixed))) {
        warning("The length of the method formals does not match the length of the parameters supplied.")
    }
    ttt.formals <- old.formals[param.free.names]
    names(ttt.formals) <- param.free.names
    new.formals <- c(ttt.formals, param.fixed)
    formals(model.req) <- new.formals
    formals(model.gen) <- new.formals
    formals(eval.priorLL) <- new.formals
    llFunc <- function(theta) {
        if (!do.call(model.req, as.list(theta))) 
            return(myRejectionLL)
        model = do.call(model.gen, as.list(theta))
        thetaLL = do.call(eval.priorLL, as.list(theta))
        result <- eval.clineLL(model) + thetaLL
        if (identical(is.finite(result), TRUE)) 
            return(result)
        return(myRejectionLL)
    }
    return(llFunc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
