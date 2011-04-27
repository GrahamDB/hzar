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
simpleModel<-makeSimpleCline1D(obs.ada.A,"free");
tailedModel<-makeTailedCline1D(obs.ada.A,"free");
print("starting repeated cline fit.");
#all.fitted.clines.d2<-list();
#iter<-1;
require(foreach)

foreach(dataset =list(obs.ak2.A,obs.pgm2.C,obs.gsr.B,obs.L5.B,obs.ada.A)) %dopar% {
  fitted.cline<-NULL;
  print("fitting:");
  print(dataset);
  try(fitted.cline<-fitClineModel(makeTailedCline1D(dataset,"none"),dataset,mcmc=1e5));
  ## if(is.null(fitted.cline)){
##     list(FALSE);
##   } else {
    fitted.cline;
##   }
  ## all.fitted.clines.d2[[iter]]<<-fitted.cline;
##   iter<<-iter+1;

} -> all.fitted.clines.d3;

  plot(cumsum(exp(sort(model.LL)))/sum(exp(sort(model.LL)))~sort(model.LL-max(model.LL)),xlim=c(-8,0),
       log="y",data=all.fitted.clines.d3[[5]]$allClines);
  lines(x=c(-4,0),y=c(0.05,0.05));
  lines(x=c(-6,-2),y=c(0.005,0.005));
  lines(x=c(-2,-2),y=c(0.01,1));
  lines(x=c(-4,-4),y=c(0.001,0.1));
