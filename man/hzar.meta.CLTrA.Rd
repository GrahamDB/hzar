\name{hzar.meta.CLTrA}
\alias{hzar.meta.CLTrA}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(hzar.meta.CLTrA)}
\format{
  The format is:
List of 4
 $ req           :function (center, width, xMin, xMax, deltaR, tauR)  
 $ prior         :function (center, width, xMin, xMax, deltaR, tauR)  
 $ func          :function (center, width, xMin, xMax, deltaR, tauR)  
 $ parameterTypes:List of 6
  ..$ center:List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "center"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ width :List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "width"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 1e-08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ xMin  :List of 2
  .. ..$ val: num 0
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "xMin"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ xMax  :List of 2
  .. ..$ val: num 10
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "xMax"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num -1e+08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ deltaR:List of 2
  .. ..$ val: num 1
  .. ..$ w  : num 1.5
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "deltaR"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 1e-08
  .. ..- attr(*, "limit.upper")= num 1e+08
  .. ..- attr(*, "realBTWN01")= logi FALSE
  ..$ tauR  :List of 2
  .. ..$ val: num 0.5
  .. ..$ w  : num 1.1
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "tauR"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 0
  .. ..- attr(*, "limit.upper")= num 1
  .. ..- attr(*, "realBTWN01")= logi TRUE
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(hzar.meta.CLTrA)
## maybe str(hzar.meta.CLTrA) ; plot(hzar.meta.CLTrA) ...
}
\keyword{datasets}
