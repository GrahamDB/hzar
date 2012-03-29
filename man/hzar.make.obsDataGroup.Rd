\name{hzar.make.obsDataGroup}
\alias{hzar.make.obsDataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Collect optimzer output from multiple models for analysis.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.make.obsDataGroup(dataGroups, obsDataGroup = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroups}{
%%     ~~Describe \code{dataGroups} here~~
}
  \item{obsDataGroup}{
%%     ~~Describe \code{obsDataGroup} here~~
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
% function (dataGroups, obsDataGroup = NULL) 
% {
%     if (identical(is.null(obsDataGroup), TRUE)) {
%         if (inherits(dataGroups, c("hzar.obsDataGroup"))) {
%             return(dataGroups)
%         }
%         if (inherits(dataGroups, c("hzar.fitRequest", "hzar.dataGroup"))) {
%             dataGroups <- list(hzar.fit2DataGroup(dataGroups))
%             obsData <- hzar.extract.obsData(dataGroups[[1]])
%             obsDataGroup <- list(data.groups = dataGroups, obsData = obsData)
%             class(obsDataGroup) <- "hzar.obsDataGroup"
%             return(obsDataGroup)
%         }
%         if (inherits(dataGroups, c("hzar.obsData"))) {
%             obsDataGroup <- list(data.groups = list(), obsData = dataGroups)
%             class(obsDataGroup) <- "hzar.obsDataGroup"
%             return(obsDataGroup)
%         }
%         if (is.list(dataGroups)) {
%             print(length(dataGroups))
%             if (inherits(dataGroups[[1]], c("hzar.fitRequest", 
%                 "hzar.dataGroup"))) {
%                 obsDataGroup <- hzar.make.obsDataGroup(dataGroups[[1]])
%             }
%             else if (inherits(dataGroups[[1]], c("hzar.obsDataGroup"))) {
%                 obsDataGroup <- dataGroups[[1]]
%             }
%             else {
%                 return(hzar.make.obsDataGroup(lapply(dataGroups, 
%                   hzar.make.obsDataGroup)))
%             }
%             if (length(dataGroups) == 1) 
%                 return(obsDataGroup)
%             otherDataGroups <- lapply(2:length(dataGroups), function(x, 
%                 y) y[[x]], y = dataGroups)
%             return(hzar.make.obsDataGroup(otherDataGroups, obsDataGroup))
%         }
%     }
  %   else if (inherits(obsDataGroup, c("hzar.obsDataGroup"))) {
  %       if (inherits(dataGroups, c("hzar.obsDataGroup"))) {
  %           return(hzar.make.obsDataGroup(dataGroups$data.groups, 
  %               obsDataGroup))
  %       }
  %       if (inherits(dataGroups, c("hzar.fitRequest", "hzar.dataGroup"))) {
  %           dataGroups <- list(hzar.fit2DataGroup(dataGroups))
  %       }
  %       if (!hzar.sameObsData(obsDataGroup, dataGroups[[1]])) 
  %           stop("All dataGroups must be from the same observation data.")
  %       if (!inherits(dataGroups[[1]], c("hzar.fitRequest", "hzar.dataGroup"))) {
  %           if (inherits(dataGroups[[1]], c("hzar.obsDataGroup"))) 
  %               return(hzar.make.obsDataGroup(hzar.make.obsDataGroup(dataGroups), 
  %                 obsDataGroup))
  %           return(hzar.make.obsDataGroup(lapply(dataGroups, 
  %               hzar.make.obsDataGroup), obsDataGroup))
  %       }
  %       modelsKnown <- as.logical(lapply(obsDataGroup$data.groups, 
  %           hzar.sameModel, fitB = dataGroups[[1]]))
  %       if (sum(modelsKnown) == 0) {
  %           obsDataGroup$data.groups <- c(obsDataGroup$data.groups, 
  %               list(dataGroups[[1]]))
  %       }
  %       else if (sum(modelsKnown) == 1) {
  %           index = which(modelsKnown)
  %           obsDataGroup$data.groups[[index]] <- hzar.dataGroup.add(obsDataGroup$data.groups[[index]], 
  %               dataGroups[[1]])
  %       }
  %       else {
  %           stop("hzar.make.obsDataGroup found a fatal error in the obsDataGroup structure.")
  %       }
  %       if (length(dataGroups) == 1) 
  %           return(obsDataGroup)
  %       otherDataGroups <- lapply(2:length(dataGroups), function(x, 
  %           y) y[[x]], y = dataGroups)
  %       return(hzar.make.obsDataGroup(otherDataGroups, obsDataGroup))
  %   }
  %   stop("Argument types unhandled by hzar.make.obsDataGroup")
  % }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
