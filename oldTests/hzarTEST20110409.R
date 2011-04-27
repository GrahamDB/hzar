print("file start");
source("hzarClasses.R");
print("classes loaded.");
mannequin.data<-read.csv("dataFiles/molecular-data-robb2.csv");
print("mannequin data loaded.");

attach(mannequin.data);

obs.ada.A<-doMolecularData1DPops(distance,p.ada.A,ada.N);
obs.ak2.A<-doMolecularData1DPops(distance,p.ak2.A,ak2.N);
obs.pgm2.C<-doMolecularData1DPops(distance,pgm2.c,pgm2.n);
obs.gsr.B<-doMolecularData1DPops(distance,p.gsr.B,gsr.N);
obs.L5.B<-doMolecularData1DPops(distance,L5.b,L5.n);

detach();
print("obs.* loaded.");
## cline.unscale <- function(cline.data, pMin,pMax){

##   cline.data$obsFreq <- (cline.data$obsFreq - pMin)/(pMax-pMin);
##   return(cline.data);
## }

## cline.logit <- function(cline.data){
##   attach(cline.data);
##   cline.data$obsFreq <- log( obsFreq / ( 1- obsFreq));
##   detach();
##   return(cline.data);
## }

## cline.samp.dist <- function(cline.data) {
##   dist.pool <- unique(sort(as.numeric(cline.data$dist)));
##   dist.init= quantile(dist.pool, probs=c(0.25,0.75));
##   potLeftOut=length(dist.pool[dist.pool <= dist.init[[1]]]);
##   potRightOut=length(dist.pool[dist.pool >= dist.init[[2]]])-1;
##   potInner=length(dist.pool[(dist.pool < dist.init[[2]])&(dist.pool > dist.init[[1]])]);
## ##   print(dist.init);
## ##   print(potLeftOut);
## ##   print(potRightOut);
## ##   print(potInner);
##   dist.len<-length(dist.pool);
##   distLeftSeries<-numeric(1);
##   distLeftSeries[1]<-dist.init[[1]];
##   distRightSeries<-numeric(1);
  
##   distRightSeries[1]<-dist.init[[2]];
##   if(potLeftOut>2)  distLeftSeries  <- c(dist.pool[2:potLeftOut],
##                                          distLeftSeries);
##   if(potRightOut>1) distRightSeries <- c(distRightSeries        ,
##                                          dist.pool[(dist.len-potRightOut):(dist.len-1)]);
##   if(potInner > 1) {
##     distLeftSeries  <- c(distLeftSeries,
##                          dist.pool[potLeftOut+1]);
##     distRightSeries <- c(dist.pool[dist.len-potRightOut-1],
##                          distRightSeries);
##   }
##   return(list(unique(distLeftSeries),unique(distRightSeries)));
## }

## cssp<-function(cline.data,do.tail="none"){
##   do.left<-FALSE;
##   do.right<-FALSE;
##   do.mirror<-FALSE;
 
##   keyValues=c("center","width","pMin","pMax")
    
##   attach(cline.data);
##   pOQ<-quantile(obsFreq,probs=c(0,0.1,0.9,1));
  
##   detach(cline.data);
##   pMinS<-as.numeric(unique(c(0,pOQ[[1]]/2,pOQ[[1]],pOQ[[2]])));
##   pMaxS<-as.numeric(unique(c(1,(1+pOQ[[4]])/2,pOQ[[4]],pOQ[[3]])));
##   bothDist<-cline.samp.dist(cline.data);
  
##   if(identical(tolower(do.tail),"mirror")){
##     do.left<-TRUE;
##     do.right<-TRUE;
##     do.mirror<-TRUE;
##     result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]],dist.right=bothDist[[2]]);
##     keyValues=c(keyValues,"deltaL","tauL");
##   } else if (identical(tolower(do.tail),"both")) {
##     do.left<-TRUE;
##     do.right<-TRUE;
##     result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]],dist.right=bothDist[[2]]);
##     keyValues=c(keyValues,"deltaL","tauL","deltaR","tauR");
##   } else if (identical(tolower(do.tail),"left")) {
##     do.left<-TRUE;
##     result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.left=bothDist[[1]]);
##     keyValues=c(keyValues,"deltaL","tauL");
##   } else if (identical(tolower(do.tail),"right")) {
##     do.right<-TRUE;
##     result<-expand.grid(pMin=pMinS,pMax=pMaxS,dist.right=bothDist[[2]]);
##     keyValues=c(keyValues,"deltaR","tauR");
##   } else {
##     result<-expand.grid(pMin=pMinS,pMax=pMaxS)
##   }
  
##   nThetas<-dim(result)[[1]];
##   result$center<-numeric(nThetas);
##   result$width <-numeric(nThetas);
##   result$deltaL <-numeric(nThetas);
##   result$tauL <-numeric(nThetas);
##   result$deltaR <-numeric(nThetas);
##   result$tauR <-numeric(nThetas);
##   for(iter in 1:nThetas){
##     tPMin<- result$pMin[[iter]];
##     tPMax<- result$pMax[[iter]];
##     myClineData=subset(cline.data,
##       (cline.data$obsFreq >tPMin)&
##       (cline.data$obsFreq <tPMax));
##     myLogitData<-cline.logit(cline.unscale(myClineData,tPMin,tPMax));
##     if(do.left){
##       myLeftData=subset(myLogitData,myLogitData$dist <= result$dist.left[[iter]]);
##       myLogitData=subset(myLogitData,myLogitData$dist >= result$dist.left[[iter]]);
##     }
##     if(do.right){
##       myRightData=subset(myLogitData,myLogitData$dist >= result$dist.right[[iter]]);
##       myLogitData=subset(myLogitData,myLogitData$dist <= result$dist.right[[iter]]);
##     }
##     lmEst<-lm(obsFreq~dist,data=myLogitData);
##     lambda=lmEst$coefficients[[2]]
##     result$width[[iter]]<- abs(0.25/lambda);
##     result$center[[iter]]<- -lmEst$coefficients[[1]]/lmEst$coefficients[[2]];
##     if(do.mirror){
##       result$center[[iter]]<- ( result$dist.left[[iter]]+ result$dist.right[[iter]])/2;
##     }
##     if(do.left){
##       lmEst<-lm(obsFreq~dist,data=myLeftData);
##       result$tauL[[iter]]<- lmEst$coefficients[[2]]/lambda;
##       result$deltaL[[iter]]<- result$center[[iter]]- result$dist.left[[iter]];
##     }

##     if(do.right){
##       lmEst<-lm(obsFreq~dist,data=myRightData);
##       result$tauR[[iter]]<- lmEst$coefficients[[2]]/lambda;
##       result$deltaR[[iter]]<-result$dist.right[[iter]]- result$center[[iter]]
##     }
##     if(do.mirror){
##       tau.mirror=( result$tauL[[iter]]+ result$tauR[[iter]])/2;
##       result$tauR[[iter]]<-tau.mirror;
##       result$tauL[[iter]]<-tau.mirror;
##     }
    
##   }
  
##   result=result[,keyValues]
##   result=na.omit(result);
##   return(result);
## } 

simpleModel<-makeSimpleCline1D(obs.ada.A,"free");
tailedModel<-makeTailedCline1D(obs.ada.A,"free");
print("starting repeated cline fit.");
all.fitted.clines.d2<-list();
iter<-1;
for(dataset in list(obs.ada.A,obs.ak2.A,obs.pgm2.C,obs.gsr.B,obs.L5.B)) {
  fitted.cline<-NULL;
  print("fitting:");
  print(dataset);
  try(fitted.cline<-fitClineModel(makeTailedCline1D(dataset,"free"),dataset,mcmc=1e6));
  if(is.null(fitted.cline)) next;
  all.fitted.clines.d2[[iter]]<<-fitted.cline;
  iter<<-iter+1;
  plot(cumsum(exp(sort(model.LL)))/sum(exp(sort(model.LL)))~sort(model.LL-max(model.LL)),xlim=c(-8,0),
       log="y",data=fitted.cline$allClines);
  lines(x=c(-4,0),y=c(0.05,0.05));
  lines(x=c(-6,-2),y=c(0.005,0.005));
  lines(x=c(-2,-2),y=c(0.01,1));
  lines(x=c(-4,-4),y=c(0.001,0.1));
}
