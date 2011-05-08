

source("hzarClasses.R");

                                            
railData<-read.csv("RailGenData.csv",header=TRUE)

## railData

obs.p.139K <- doMolecularData1DPops(railData$dist, railData$p.139K, railData$p.139.N)

## obs.p.139K


model3.cline.p.139K <- makeCline1D(obs.p.139K, scaling="free", tail="right")

model3.fit.cline.p.139K <- fitClineModel(model3.cline.p.139K, obs.p.139K, mcmc=1e6)

## error:

##               center        width     deltaR         tauR model.LL
## center    0.782549718 -0.009108033  3.2003879 -0.026109313      NaN
## width    -0.009108033  0.020368379  0.9663805 -0.001178561      NaN
## deltaR    3.200387906  0.966380524 83.8164639 -0.387211881      NaN
## tauR     -0.026109313 -0.001178561 -0.3872119  0.003277250      NaN
## model.LL          NaN          NaN        NaN          NaN      NaN
## [1] "C"
##       center width deltaR tauR
## center    NaN   NaN    NaN  NaN
## width     NaN   NaN    NaN  NaN
## deltaR    NaN   NaN    NaN  NaN
## tauR      NaN   NaN    NaN  NaN
## Error in if (!do.call(meta.model$req, as.list(theta))) return(myRejectionLL) :
##  missing value where TRUE/FALSE needed



## something odd                          
##      center     width       pMin   pMax    deltaR         tauR   model.LL
## 51 6.053529 0.2574508 0.07894737 1.0000  2.546471 8.226173e-02       -Inf
## 52 6.294688 0.2379888 0.12868421 1.0000  2.305312 7.608287e-02       -Inf
## 53 6.154532 0.2504819 0.10381579 1.0000  2.445468 8.005530e-02       -Inf
## 54 6.383559 0.2282365 0.13868421 1.0000  2.216441 7.297337e-02       -Inf
