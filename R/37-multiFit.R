
hzar.doFit.multi <- function(mFitR,doPar=FALSE,inOrder=TRUE){
  fitR=NULL;
  res <- list();
  if(doPar){
    res <- foreach(fitR=mFitR,.inorder=inOrder) %dopar% {
      hzar.doFit(fitR);
    }
  }else {
    res <- foreach(fitR=mFitR,.inorder=inOrder) %do% {
      hzar.doFit(fitR);
    }
  }
  res
}

hzar.doChain.multi <- function(mFitR,doPar=FALSE,inOrder=TRUE,...){
  fitR=NULL;
  res <- list();
  if(doPar){
    res <- foreach(fitR=mFitR,.inorder=inOrder) %dopar% {
      hzar.chain.doSeq(fitR,...);
    }
  }else {
    res <- foreach(fitR=mFitR,.inorder=inOrder) %do% {
      hzar.chain.doSeq(fitR,...);
    }
  }
  res
}
