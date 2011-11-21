
model.addReqClause <- function(meta.model,newClause){
  oldClause <- body(meta.model$req)[[2]][[2]];
  body(meta.model$req)[[2]][[2]]<-substitute(clause1 & clause2,list(clause1=oldClause,clause2=newClause));
  return(meta.model);
}

hzar.model.addMaxWidth <- function(meta.model,maxWidth){
  clause <- substitute(width<V,list(V=maxWidth));
  attr(meta.model$parameterTypes$width,"limit.upper") <- maxWidth;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addMaxCenter <- function(meta.model,maxCenter){
  clause <- substitute(center<V,list(V=maxCenter));
  attr(meta.model$parameterTypes$center,"limit.upper") <- maxCenter;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addMaxDelta <- function(meta.model,maxDelta){
  if(identical(attr(meta.model,"tails"),"none"))
    return(meta.model);
  
  clauseL <- substitute(deltaL<V,list(V=maxDelta));
  clauseM <- substitute(deltaM<V,list(V=maxDelta));
  clauseR <- substitute(deltaR<V,list(V=maxDelta));
  assignLimit <- function(param){
    attr(meta.model$parameterTypes[[param]],"limit.upper") <<- maxDelta;
  }
  if(identical(attr(meta.model,"tails"),"both")){
    assignLimit("deltaL");
    assignLimit("deltaR");
    return(model.addReqClause(model.addReqClause(meta.model,clauseL),clauseR));
  }
  if(identical(attr(meta.model,"tails"),"right")){
    assignLimit("deltaR");
    return(model.addReqClause(meta.model,clauseR));
  }
  if(identical(attr(meta.model,"tails"),"left")){
    assignLimit("deltaL");
    return(model.addReqClause(meta.model,clauseL));
  }
  
  if(identical(attr(meta.model,"tails"),"mirror")){
    assignLimit("deltaM");
    return(model.addReqClause(meta.model,clauseM));
  }
  
}
hzar.model.addMinCenter <- function(meta.model,minCenter){
  clause <- substitute(center>V,list(V=minCenter));
  attr(meta.model$parameterTypes$center,"limit.lower") <- minCenter;
  return(model.addReqClause(meta.model,clause));
}

hzar.model.addCenterRange <- function(meta.model,low,high){
  return(hzar.model.addMaxCenter(hzar.model.addMinCenter(meta.model,low),high))
}

hzar.model.addBoxReq <- function(meta.model,low,high){
  return(hzar.model.addCenterRange(hzar.model.addMaxDelta(hzar.model.addMaxWidth(meta.model,high-low),high-low),low,high));
}
