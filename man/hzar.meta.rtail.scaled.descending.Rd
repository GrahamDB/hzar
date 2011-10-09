\name{hzar.meta.rtail.scaled.descending}
\alias{hzar.meta.rtail.scaled.descending}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(hzar.meta.rtail.scaled.descending)}
\format{
  The format is:
List of 4
 $ req           :function (center, width, pMin, pMax, deltaR, tauR)  
 $ prior         :function (center, width, pMin, pMax, deltaR, tauR)  
 $ func          :function (center, width, pMin, pMax, deltaR, tauR)  
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
  ..$ pMin  :List of 2
  .. ..$ val: num 0
  .. ..$ w  : num 1.1
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "pMin"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 0
  .. ..- attr(*, "limit.upper")= num 1
  .. ..- attr(*, "realBTWN01")= logi TRUE
  ..$ pMax  :List of 2
  .. ..$ val: num 1
  .. ..$ w  : num 1.1
  .. ..- attr(*, "class")= chr "clineParameter"
  .. ..- attr(*, "param")= chr "pMax"
  .. ..- attr(*, "fixed")= logi FALSE
  .. ..- attr(*, "limit.lower")= num 0
  .. ..- attr(*, "limit.upper")= num 1
  .. ..- attr(*, "realBTWN01")= logi TRUE
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
 - attr(*, "class")= chr "clineMetaModel"
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
data(hzar.meta.rtail.scaled.descending)
## maybe str(hzar.meta.rtail.scaled.descending) ; plot(hzar.meta.rtail.scaled.descending) ...
}
\keyword{datasets}
