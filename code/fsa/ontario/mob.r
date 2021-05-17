source('../config.r')
source('data.r')
source('plot.r')

do.mobility = function(do='plots'){
  pop = load.fsa.pop(age=FALSE)
  pop$decile = as.factor(pop$decile)
  X = load.fsa.mob()
  # aggregate over visited
  S = merge(
    aggregate(devices.home~FSA+month,X,mean,na.rm=TRUE,drop=FALSE),
    aggregate(devices.travel~FSA+month,X,sum,na.rm=TRUE,drop=FALSE))
  S = merge(S,merge(pop,load.fsa.smartphones()),all.x=TRUE)
  S = S[order(S$month,S$FSA),] # reshape for implicit repeat along FSA
  # get reference numbers of devices per FSA
  idx.ref = S$month %in% mo.ref
  dt.fsa = aggregate(devices.total~FSA+month,S,sum,na.rm=TRUE,drop=FALSE)$devices.total
  dt.g   = aggregate(devices.total~decile+month,S,sum,na.rm=TRUE,drop=FALSE)$devices.total
  dh.ref = aggregate(devices.home~FSA,S[idx.ref,],mean,drop=FALSE)$devices.home
  dt.ref = aggregate(devices.travel~FSA,S[idx.ref,],mean,drop=FALSE)$devices.travel
  X = X[order(X$month,X$FSA.visited,X$FSA),] # reshape for implicit repeat along FSA
  X$devices.total  = S$devices.total
  X$devices.travel = S$devices.trave
  # TODO: better dealing with missing propoprtion (?)
  X$travel.total = X$devices.travel * (1 + OR.travel.unobs * (X$devices.total - X$devices.home) / dh.ref) # TEMP
  # X$travel.total = X$devices.travel * (1 + OR.travel.unobs * (X$devices.total / dh.ref - 1)) # TEMP
  # X$travel.total = X$devices.travel * (1 + OR.travel.unobs * (X$devices.total / X$devices.home - 1)) # TEMP
  G.h = pop[,c('FSA','decile')]; G.v = G.h; colnames(G.v) = c('FSA.visited','decile.visited')
  if (do == 'plots'){
    S$devices.ratio = S$devices.travel / S$devices.home
    S$travel.red = S$devices.travel / dt.ref
    S$home.red   = S$devices.home   / dh.ref
    S$ratio.red  = S$travel.red     / S$home.red
    S$travel.rx  = S$devices.travel / dh.ref
    print(summary(S[!idx.ref,]$ratio.red))
    print(summary(100*S[ idx.ref,]$devices.home/S[ idx.ref,]$devices.total))
    print(summary(100*S[!idx.ref,]$devices.home/S[!idx.ref,]$devices.total))
    plot.device.density(S,'devices.home',  xmax=1500) + labs(y='Month',x='Devices at Home'); ggsave(figname('devices-month-home'),width=5,heigh=5)
    plot.device.density(S,'devices.travel',xmax=1500) + labs(y='Month',x='Device Visits');   ggsave(figname('devices-month-travel'),width=5,heigh=5)
    plot.device.density(S,'devices.ratio', xmax=2)    + labs(y='Month',x='Visits / Device'); ggsave(figname('devices-month-ratio'),width=5,heigh=5)
    plot.device.density(S,'home.red',  xmax=2) + labs(y='Month',x='Home Reduction');         ggsave(figname('reduction-month-home'),width=5,heigh=5)
    plot.device.density(S,'travel.red',xmax=2) + labs(y='Month',x='Travel Reduction');       ggsave(figname('reduction-month-travel'),width=5,heigh=5)
    plot.device.density(S,'ratio.red', xmax=2) + labs(y='Month',x='Travel Reduction / Home Reduction'); ggsave(figname('reduction-month-ratio'),width=5,heigh=5)
    plot.device.density(S,'travel.rx', xmax=2) + labs(y='Month',x='Travel / Reference Home'); ggsave(figname('reduction-month-travel-x'),width=5,heigh=5)
    X.h = aggregate(travel.total~FSA+month,X,sum,na.rm=TRUE,drop=FALSE)
    X.h$travel.prop = X.h$travel.total / dt.fsa
    X.h = merge(X.h,G.h,all.x=TRUE)
    plot.device.density(X.h,'travel.prop',y='month', xmax=1.2); ggsave(figname('prop-travel-month'), width=5,height=5)
    plot.device.density(X.h,'travel.prop',y='decile',xmax=1.2); ggsave(figname('prop-travel-decile'),  width=5,height=5)
  }
  if (do == 'csv'){
    X = merge(merge(G,G.h,all.x=TRUE),G.v,all.x=TRUE) # map deciles TODO: speed up by reshape?
    X.gg = aggregate(travel.total~decile+decile.visited+month,X,sum,na.rm=TRUE,drop=FALSE)
    X.gg$travel.prop = X.gg$travel.total / dt.g
    X.gg$travel.total = NULL
    X.gg = X.gg[order(X.gg$month,X.gg$decile,X.gg$decile.visited),]
    write.csv(X.gg,root.path('data','fsa','mobility.csv'),row.names=FALSE)
  }
}

do.mobility()