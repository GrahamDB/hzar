\name{cssp}
\alias{cssp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cssp(cline.data, do.tail = "none", pMinS = NULL, pMaxS = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cline.data}{
%%     ~~Describe \code{cline.data} here~~
}
  \item{do.tail}{
%%     ~~Describe \code{do.tail} here~~
}
  \item{pMinS}{
%%     ~~Describe \code{pMinS} here~~
}
  \item{pMaxS}{
%%     ~~Describe \code{pMaxS} here~~
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
function (cline.data, do.tail = "none", pMinS = NULL, pMaxS = NULL) 
{
    do.left <- FALSE
    do.right <- FALSE
    do.mirror <- FALSE
    keyValues = c("center", "width", "pMin", "pMax")
    attach(cline.data)
    pOQ <- quantile(unique(sort(obsFreq)), probs = c(0, 0.1, 
        0.9, 1))
    detach(cline.data)
    if (is.null(pMinS)) {
        pMinS <- as.numeric(unique(c(0, pOQ[[1]]/2, pOQ[[1]], 
            pOQ[[2]])))
        if (length(pMaxS) < 3) {
            pMinS <- as.numeric(unique(c(pMinS, (pOQ[[1]] + pOQ[[2]])/2, 
                0.01 + pOQ[[2]])))
        }
    }
    else {
        pMinS <- as.numeric(pMinS)
        if (length(pMinS) == 0) {
            warning("Malformed pMinS, using 0 instead")
            pMinS <- as.numeric(c(0))
        }
    }
    if (is.null(pMaxS)) {
        pMaxS <- as.numeric(unique(c(1, (1 + pOQ[[4]])/2, pOQ[[4]], 
            pOQ[[3]])))
        if (length(pMaxS) < 3) {
            pMaxS <- as.numeric(unique(c(pMaxS, (pOQ[[3]] + pOQ[[4]])/2, 
                pOQ[[3]] - 0.01)))
        }
    }
    else {
        pMaxS <- as.numeric(pMaxS)
        if (length(pMaxS) == 0) {
            warning("Malformed pMaxS, using 1 instead")
            pMaxS <- as.numeric(c(1))
        }
    }
    bothDist <- cline.samp.dist(cline.data)
    if (identical(tolower(do.tail), "mirror")) {
        do.left <- TRUE
        do.right <- TRUE
        do.mirror <- TRUE
        result <- expand.grid(pMin = pMinS, pMax = pMaxS, dist.left = bothDist[[1]], 
            dist.right = bothDist[[2]])
        keyValues = c(keyValues, "deltaM", "tauM")
    }
    else if (identical(tolower(do.tail), "both")) {
        do.left <- TRUE
        do.right <- TRUE
        result <- expand.grid(pMin = pMinS, pMax = pMaxS, dist.left = bothDist[[1]], 
            dist.right = bothDist[[2]])
        keyValues = c(keyValues, "deltaL", "tauL", "deltaR", 
            "tauR")
    }
    else if (identical(tolower(do.tail), "left")) {
        do.left <- TRUE
        result <- expand.grid(pMin = pMinS, pMax = pMaxS, dist.left = bothDist[[1]])
        keyValues = c(keyValues, "deltaL", "tauL")
    }
    else if (identical(tolower(do.tail), "right")) {
        do.right <- TRUE
        result <- expand.grid(pMin = pMinS, pMax = pMaxS, dist.right = bothDist[[2]])
        keyValues = c(keyValues, "deltaR", "tauR")
    }
    else {
        result <- expand.grid(pMin = pMinS, pMax = pMaxS)
    }
    nThetas <- dim(result)[[1]]
    result$center <- numeric(nThetas)
    result$width <- numeric(nThetas)
    result$deltaL <- numeric(nThetas)
    result$tauL <- numeric(nThetas)
    result$deltaR <- numeric(nThetas)
    result$tauR <- numeric(nThetas)
    result$deltaM <- numeric(nThetas)
    result$tauM <- numeric(nThetas)
    for (iter in 1:nThetas) {
        tPMin <- result$pMin[[iter]]
        tPMax <- result$pMax[[iter]]
        myClineData <- cline.data
        myClineData$obsFreq[myClineData$obsFreq < (1e-04 + tPMin)] <- rep(1e-04 + 
            tPMin, length(which(myClineData$obsFreq < (1e-04 + 
            tPMin))))
        myClineData$obsFreq[myClineData$obsFreq > (-1e-04 + tPMax)] <- rep(tPMax - 
            1e-04, length(which(myClineData$obsFreq > (tPMax - 
            1e-04))))
        myLogitData <- cline.logit(cline.unscale(myClineData, 
            tPMin, tPMax))
        myFullLogitData <- myLogitData
        if (do.left) {
            myLeftData = subset(myLogitData, myLogitData$dist <= 
                result$dist.left[[iter]])
            myLogitData = subset(myLogitData, myLogitData$dist >= 
                result$dist.left[[iter]])
        }
        if (do.right) {
            myRightData = subset(myLogitData, myLogitData$dist >= 
                result$dist.right[[iter]])
            myLogitData = subset(myLogitData, myLogitData$dist <= 
                result$dist.right[[iter]])
        }
        lmEst <- NULL
        try(lmEst <- lm(obsFreq ~ dist, data = myLogitData))
        if (is.null(lmEst)) {
            try(lmEst <- lm(obsFreq ~ dist, data = myFullLogitData))
            if (is.null(lmEst)) {
                result$width[[iter]] <- NA
                next
            }
        }
        lambda = lmEst$coefficients[[2]]
        result$width[[iter]] <- abs(0.25/lambda)
        result$center[[iter]] <- -lmEst$coefficients[[1]]/lmEst$coefficients[[2]]
        if (do.mirror) {
            result$center[[iter]] <- (result$dist.left[[iter]] + 
                result$dist.right[[iter]])/2
        }
        if (do.left) {
            if (result$center[[iter]] < result$dist.left[[iter]]) 
                result$center[[iter]] = result$dist.left[[iter]]
            lmEst <- NULL
            try(lmEst <- lm(obsFreq ~ dist, data = myLeftData))
            if (is.null(lmEst)) {
                result$tauL[[iter]] <- 0
            }
            else {
                propose.tau <- lmEst$coefficients[[2]]/lambda
                propose.tau = ifelse(propose.tau < 0, 1e-04 + 
                  propose.tau * 1e-05, propose.tau)
                propose.tau = ifelse(propose.tau >= 1, 1 + (propose.tau - 
                  1) * 1e-05 - 1e-04, propose.tau)
                propose.tau = ifelse(propose.tau < 0, 0, propose.tau)
                propose.tau = ifelse(propose.tau > 1, 1, propose.tau)
                result$tauL[[iter]] <- propose.tau
            }
            result$deltaL[[iter]] <- result$center[[iter]] - 
                result$dist.left[[iter]]
        }
        if (do.right) {
            if (result$center[[iter]] > result$dist.right[[iter]]) 
                result$center[[iter]] = result$dist.right[[iter]]
            lmEst <- NULL
            try(lmEst <- lm(obsFreq ~ dist, data = myRightData))
            if (is.null(lmEst)) {
                result$tauR[[iter]] <- 0
            }
            else {
                propose.tau <- lmEst$coefficients[[2]]/lambda
                propose.tau = ifelse(propose.tau < 0, 1e-04 + 
                  propose.tau * 1e-05, propose.tau)
                propose.tau = ifelse(propose.tau >= 1, 1 + (propose.tau - 
                  1) * 1e-05 - 1e-04, propose.tau)
                propose.tau = ifelse(propose.tau < 0, 0, propose.tau)
                propose.tau = ifelse(propose.tau > 1, 1, propose.tau)
                result$tauR[[iter]] <- propose.tau
            }
            result$deltaR[[iter]] <- result$dist.right[[iter]] - 
                result$center[[iter]]
        }
        if (do.mirror) {
            tau.mirror = (result$tauL[[iter]] + result$tauR[[iter]])/2
            result$tauR[[iter]] <- tau.mirror
            result$tauL[[iter]] <- tau.mirror
            result$tauM[[iter]] <- tau.mirror
            result$deltaM[[iter]] <- result$deltaL[[iter]]
        }
    }
    result = result[, keyValues]
    result = na.omit(result)
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
