source('data.r')
source('plot.r')

# Objective: 1. estimate mobility matrix (decile) from cell-phone count data (FSA)
#            2. plot some marginal distributions of the mobility data

g.merge = function(X,pop,case){
  cols = c('FSA','decile')
  g = pop[,cols]
  colnames(g) = switch(case,'h'=cols,'v'=paste0(cols,'.visited'))
  return(merge(X,g,all.x=TRUE))
}

plot.mobility.margins = function(X,S,pop,idx.ref,f='mobility'){
  dh.ref = aggregate(devices.home ~FSA,S[idx.ref,],mean,drop=FALSE)$devices.home
  dv.ref = aggregate(devices.visit~FSA,S[idx.ref,],mean,drop=FALSE)$devices.visit
  S$devices.vt.ht = S$devices.visit / S$devices.home
  S$devices.vt.v0 = S$devices.visit / dv.ref
  S$devices.ht.h0 = S$devices.home  / dh.ref
  S$devices.vr.hr = S$devices.vt.v0 / S$devices.ht.h0
  S$devices.vt.h0 = S$devices.visit / dh.ref
  print(summary(S[!idx.ref,]$devices.vr.hr))
  print(summary(100*S[ idx.ref,]$devices.home / S[ idx.ref,]$devices.total))
  print(summary(100*S[!idx.ref,]$devices.home / S[!idx.ref,]$devices.total))
  print(summary(100*(S$pop - S$devices.total) / S$pop))
  g = plot.device.density(list(
    'Home Reduction'         = data.frame(month=S$month,Ratio=S$devices.ht.h0),
    'Visit Reduction'        = data.frame(month=S$month,Ratio=S$devices.vt.v0),
    'Visit / Reference Home' = data.frame(month=S$month,Ratio=S$devices.vt.h0)
  ),'Ratio',xmax=2); ggsave(figname('devices-panel','mx'),w=15,h=5)
  plot.device.density(S,'devices.home', xmax=1500) + xlab('Devices at Home'); ggsave(figname('devices-n-month-home', f),width=5,heigh=5)
  plot.device.density(S,'devices.visit',xmax=1500) + xlab('Device Visits');   ggsave(figname('devices-n-month-visit',f),width=5,heigh=5)
  plot.device.density(S,'devices.vt.ht',xmax=2)    + xlab('Visits / Device'); ggsave(figname('devices-n-month-ratio',f),width=5,heigh=5)
  plot.device.density(S,'devices.ht.h0',xmax=2) + xlab('Home Reduction');     ggsave(figname('devices-r-month-home', f),width=5,heigh=5)
  plot.device.density(S,'devices.vt.v0',xmax=2) + xlab('Visit Reduction');    ggsave(figname('devices-r-month-visit',f),width=5,heigh=5)
  plot.device.density(S,'devices.vr.hr',xmax=2) + xlab('Visit Reduction / Home Reduction'); ggsave(figname('devices-r-month-ratio',f),width=5,heigh=5)
  plot.device.density(S,'devices.vt.h0',xmax=2) + xlab('Visit / Reference Home'); ggsave(figname('devices-r-month-rx',f),width=5,heigh=5)
  plot.device.density(S[idx.ref,],'devices.vt.h0',y='decile',xmax=2) + xlab('Visit / Reference Home'); ggsave(figname('devices-r-decile-rx',f),width=5,heigh=5)
  X.h = aggregate(visit.total~FSA+month,X,sum,na.rm=TRUE,drop=FALSE)
  X.h = merge(X.h,pop,all.x=TRUE)
  X.h$visit.prop = X.h$visit.total / X.h$pop
  plot.device.density(X.h,'visit.prop',xmax=2); ggsave(figname('travel-p-month',f),width=5,heigh=5)
  plot.device.density(X.h,'visit.prop',y='decile',xmax=2); ggsave(figname('travel-p-decile',f),width=5,heigh=5)
  # TODO: plot monthly mobility here from mixing
}

save.mobility = function(X,S,pop){
  X = g.merge(g.merge(X,pop,'h'),pop,'v') # WARNING: this is slow
  pop.g = aggregate(pop~decile+month,S,sum,na.rm=TRUE,drop=FALSE)$pop
  X.gg = aggregate(visit.total~decile+decile.visited+month,X,sum,na.rm=TRUE,drop=FALSE)
  X.gg$visit.prop = X.gg$visit.total / pop.g
  X.gg$visit.total = NULL
  X.gg = X.gg[order(X.gg$month,X.gg$decile,X.gg$decile.visited),]
  write.csv(X.gg,root.path('data','mix','mobility_decile.csv'),row.names=FALSE)
}

main.mobility = function(do.plot=TRUE,do.csv=TRUE){
  pop = load.fsa.pop(age=FALSE)
  pop$decile = as.factor(pop$decile)
  X = load.fsa.mob()
  # aggregate over visited
  S = merge(
    aggregate(devices.home ~FSA+month,X,mean,na.rm=TRUE,drop=FALSE),
    aggregate(devices.visit~FSA+month,X,sum, na.rm=TRUE,drop=FALSE))
  S = merge(S,merge(pop,load.fsa.smartphones()),all.x=TRUE)
  S = S[order(S$month,S$FSA),] # reshape for implicit repeat along FSA
  # get reference numbers of devices per FSA
  idx.ref = S$month %in% config$t.ref
  dh.ref = aggregate(devices.home ~FSA,S[idx.ref,],mean,drop=FALSE)$devices.home
  X = X[order(X$month,X$FSA.visited,X$FSA),] # reshape for implicit repeat along FSA
  X$visit.total = X$devices.visit / dh.ref * (1 # "extrapolate" the unobserved people
    + config$phi['unobs.device'] * (S$devices.total - dh.ref)
    + config$phi['no.device'] * (S$pop - S$devices.total))
  if (do.plot){
    plot.mobility.margins(X,S,pop,idx.ref)
  }
  if (do.csv){
    save.mobility(X,S,pop)
  }
}
