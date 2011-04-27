source("hzarClasses.R");
mannequin.data<-read.csv("dataFiles/molecular-data-robb2.csv");

attach(mannequin.data);
obsAdaA<-doMolecularData1DPops(distance,p.ada.A,ada.N);
detach();
simpleModel<-makeSimpleCline1D(obsAdaA,"free");
tailedModel<-makeTailedCline1D(obsAdaA,"free");

a.subsets<-list(c(1,2),c(1,3),c(1,4),c(1,5),c(1,6),c(1,7),c(1,8),
                c(2,3),c(2,4),c(2,5),c(2,6),c(2,7),c(2,8),
                c(3,4),c(3,5),c(3,6),c(3,7),c(3,8),
                c(4,5),c(4,6),c(4,7),c(4,8),
                c(5,6),c(5,7),c(5,8),
                c(6,7),c(6,8),
                c(7,8));
b.subsets<-list(c(3,4,5,6,7,8),c(2,4,5,6,7,8),c(2,3,5,6,7,8),
                c(2,3,4,6,7,8),c(2,3,4,5,7,8),c(2,3,4,5,6,8),c(2,3,4,5,6,7),
                c(1,4,5,6,7,8),c(1,3,5,6,7,8),c(1,3,4,6,7,8),
                c(1,3,4,5,7,8),c(1,3,4,5,6,8),c(1,3,4,5,6,7),
                c(1,2,5,6,7,8),c(1,2,4,6,7,8),
                c(1,2,4,5,7,8),c(1,2,4,5,6,8),c(1,2,4,5,6,7),
                c(1,2,3,6,7,8),c(1,2,3,5,7,8),c(1,2,3,5,6,8),c(1,2,3,5,6,7),
                c(1,2,3,4,7,8),c(1,2,3,4,6,8),c(1,2,3,4,6,7),
                c(1,2,3,4,5,8),c(1,2,3,4,5,7),c(1,2,3,4,5,6));
attach(tailedModel);
param.full<-splitParameters(parameterTypes)
results.out<-list(1:28);                
results.CC<-list(1:28); 
for(iter in 1:28){
  a.set<-a.subsets[[iter]];
  b.set<-b.subsets[[iter]];
tttFunc.full<-function(ttt) ifelse(do.call(req,as.list(c(ttt,
                                                         param.full$init[b.set],
                                                         param.full$fixed))),
                                   obsAdaA$model.LL(do.call(func,
                                                            as.list(c(ttt,
                                                                      param.full$init[b.set],
                                                                      param.full$fixed)))),
                                   -1e8);
  print(param.full$init[a.set]);
  opt.out<-NULL
  CC<-NULL
try(optim(param.full$init[a.set],tttFunc.full,
      control=list(fnscale=-1,trace=0,maxit=1e3),
      method="L-BFGS-B",upper=param.full$upper[a.set],lower=param.full$lower[a.set],
      hessian=TRUE) -> opt.out );
results.out[[iter]]<-opt.out;
try(print(opt.out$hessian));
try(CC <- chol(-1 * opt.out$hessian));
results.CC[[iter]]<-CC;
  try(print(CC));
}
detach();
