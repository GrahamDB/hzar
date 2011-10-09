\name{setupMoleCenterClineParameters}
\alias{setupMoleCenterClineParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setupMoleCenterClineParameters(myModel, scaling, x = NULL, y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myModel}{
%%     ~~Describe \code{myModel} here~~
}
  \item{scaling}{
%%     ~~Describe \code{scaling} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (myModel, scaling, x = NULL, y = NULL) 
{
    pTnames <- names(myModel$parameterTypes)
    if ("pMin" \%in\% pTnames) 
        mdlMin <- "pMin"
    if ("xMin" \%in\% pTnames) 
        mdlMin <- "xMin"
    if ("pMax" \%in\% pTnames) 
        mdlMax <- "pMax"
    if ("xMax" \%in\% pTnames) 
        mdlMax <- "xMax"
    if (scaling == "none") {
        attr(myModel$parameterTypes[[mdlMin]], "fixed") <- TRUE
        attr(myModel$parameterTypes[[mdlMax]], "fixed") <- TRUE
        myModel$parameterTypes[[mdlMin]]$val <- 0
        myModel$parameterTypes[[mdlMax]]$val <- 1
    }
    else if (scaling == "fixed") {
        attr(myModel$parameterTypes[[mdlMin]], "fixed") <- TRUE
        attr(myModel$parameterTypes[[mdlMax]], "fixed") <- TRUE
        if (!is.null(y)) {
            myModel$parameterTypes[[mdlMin]]$val <- min(y)
            myModel$parameterTypes[[mdlMax]]$val <- max(y)
        }
    }
    else if (scaling == "free") {
        attr(myModel$parameterTypes[[mdlMin]], "fixed") <- FALSE
        attr(myModel$parameterTypes[[mdlMax]], "fixed") <- FALSE
        if (!is.null(y)) {
            myModel$parameterTypes[[mdlMin]]$val <- min(y)
            myModel$parameterTypes[[mdlMax]]$val <- max(y)
            junk <- hzar.suggestionFunc1D[[mdlMin]](x, y)
            attr(myModel$parameterTypes[[mdlMin]], "limit.lower") <- junk[[1]]
            attr(myModel$parameterTypes[[mdlMin]], "limit.upper") <- junk[[2]]
            junk <- hzar.suggestionFunc1D[[mdlMax]](x, y)
            attr(myModel$parameterTypes[[mdlMax]], "limit.lower") <- junk[[1]]
            attr(myModel$parameterTypes[[mdlMax]], "limit.upper") <- junk[[2]]
        }
    }
    else {
        stop(paste("Scaling type", scaling, "unrecignized. Please use none, fixed, or free."))
    }
    if (!is.null(x)) {
        qX <- quantile(x, probs = c(0.25, 0.5, 0.75))
        myModel$parameterTypes$center$val <- qX[[2]]
        myModel$parameterTypes$width$val <- qX[[3]] - qX[[1]]
        junk <- hzar.suggestionFunc1D$center(x, y)
        attr(myModel$parameterTypes$center, "limit.lower") <- junk[[1]]
        attr(myModel$parameterTypes$center, "limit.upper") <- junk[[2]]
        junk <- hzar.suggestionFunc1D$width(x, y)
        attr(myModel$parameterTypes$width, "limit.lower") <- junk[[1]]
        attr(myModel$parameterTypes$width, "limit.upper") <- junk[[2]]
        index <- "deltaR"
        if (index \%in\% pTnames) {
            junk <- hzar.suggestionFunc1D[[index]](x, y)
            attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
            attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
        }
        index <- "deltaM"
        if (index \%in\% pTnames) {
            junk <- hzar.suggestionFunc1D[[index]](x, y)
            attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
            attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
        }
        index <- "deltaL"
        if (index \%in\% pTnames) {
            junk <- hzar.suggestionFunc1D[[index]](x, y)
            attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
            attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
        }
    }
    index <- "tauR"
    if (index \%in\% pTnames) {
        junk <- hzar.suggestionFunc1D[[index]](x, y)
        attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
        attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
    }
    index <- "tauM"
    if (index \%in\% pTnames) {
        junk <- hzar.suggestionFunc1D[[index]](x, y)
        attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
        attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
    }
    index <- "tauL"
    if (index \%in\% pTnames) {
        junk <- hzar.suggestionFunc1D[[index]](x, y)
        attr(myModel$parameterTypes[[index]], "limit.lower") <- junk[[1]]
        attr(myModel$parameterTypes[[index]], "limit.upper") <- junk[[2]]
    }
    return(myModel)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
