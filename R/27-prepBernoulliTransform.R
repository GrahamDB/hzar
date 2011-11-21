
hzar.dBernoulli.LL <- function(values,locations,getMax=FALSE,getProbs=FALSE){
  midpoints<-function(x) {
    y<-as.numeric(sort(unique(x)));
    return((y[y>min(y)]+y[y<max(y)])/2);
  }
  mVal<-midpoints(values);
  bernoulli.LL <- function(nSuccess,nTotal){
    nFailures<-nTotal-nSuccess;
    hSuccess<-ifelse(nSuccess>0,nSuccess*log(nTotal/nSuccess),0);
    hFailure<-ifelse(nFailures>0,nFailures*log(nTotal/nFailures),0);
    return(hSuccess+hFailure);
  }
  pL<-lapply(unique(locations),function(x) as.numeric(values[locations==x]));
  bLLSum<- function(junk,procList){
    nSamples=as.numeric(lapply(procList,length))
    return(sum(bernoulli.LL(as.numeric(lapply(procList,
                                              function(x,cutValue)
                                              sum(x<cutValue),
                                              cutValue=junk)),
                            nSamples)));
  }
  if(getProbs){
    junk<-as.numeric(lapply(mVal, bLLSum, pL ));
    cutValue<-mVal[junk==max(junk)][1];
    nSamples<-as.numeric(lapply(pL,length));
    nLow<-as.numeric(lapply(pL,function(x,cutValue) sum(x<cutValue), cutValue) )
    return(data.frame(locID=unique(locations),
                      pObs=nLow/nSamples,
                      nLow=nLow,
                      nSamples=nSamples,
                      cutValue=cutValue));
  }
  if(getMax){
    junk<-as.numeric(lapply(mVal, bLLSum, pL ));
    return(mVal[junk==max(junk)][1]);
  }
  return(as.numeric(lapply(mVal, bLLSum, pL )));
}

hzar.makeTraitObsData <- function(distOfLocation,locationOfValue,values){
  probs <- hzar.dBernoulli.LL(values,locationOfValue,getProbs=TRUE)
  return(hzar.doMolecularData1DPops(distance=distOfLocation[probs$locID],
                                  pObs=probs$pObs,
                                  nSamples=probs$nSamples));
}

hzar.doMorphoSets <- function(traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol){
  distOfLocation <- tDist[[tDDistCol]];
  names(distOfLocation) <- tDist[[tDLocCol]];
  
  res <- lapply(traitNames,
                function(id) {
                  tVal <- tValues[!is.na(tValues[[id]]),c(tVLocCol,id)]
                  
                  return(hzar.makeTraitObsData(distOfLocation,
                                        locationOfValue=tVal[[tVLocCol]],
                                        values=tVal[[id]]));
                  }
                )
  names(res) <- traitNames;
  return(res);
}
