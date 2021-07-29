source('data.r')
source('plot.r')

# Objective: 1. estimate mobility matrix (decile) from cell-phone count data (FSA)
#            2. plot some marginal distributions of the mobility data

merge.mobility.margin = function(X,H){
  # load some needed stuff
  phones = load.fsa.smartphones()
  pop = load.fsa.pop(age=FALSE)
  pop$decile = as.factor(pop$decile)
  # aggregate over visited
  S = merge(
    aggregate(devices.home ~FSA+month,X,mean,na.rm=TRUE,drop=FALSE),
    aggregate(devices.visit~FSA+month,X,sum, na.rm=TRUE,drop=FALSE))
  S = merge(S,merge(pop,phones),all.x=TRUE)
  S = merge(S,H,all.x=TRUE)
  return(S[order(S$month,S$FSA),]) # reshape for implicit repeat along FSA
}

gen.mobility.refs = function(S,idx.ref){
  return(list(
    devices.home  = aggregate(devices.home ~FSA,S[idx.ref,],mean,drop=FALSE)$devices.home,
    devices.visit = aggregate(devices.visit~FSA,S[idx.ref,],mean,drop=FALSE)$devices.visit,
    away.time     = aggregate(away.time    ~FSA,S[idx.ref,],mean,drop=FALSE)$away.time
  ))
}

add.within.fsa.mobility = function(X,S,SR){
  X = X[order(X$month,X$FSA.visited,X$FSA),]
  idx.self = X$FSA == X$FSA.visited
  X[idx.self,]$devices.home = S$devices.home # was NA
  X[idx.self,]$devices.visit =
    pmin(S$away.ratio, 1.5) * pmax(0, SR$devices.home - S$devices.visit)
  return(X)
}

mobility.scale.up = function(X,S,SR){
  X = X[order(X$FSA.visited,X$month,X$FSA),]
  X$visit.total = X$devices.visit / SR$devices.home * (1 # "extrapolate" the unobserved people
    + config$phi['unobs.device'] * (S$devices.total - SR$devices.home)
    + config$phi['no.device'] * (S$pop - S$devices.total))
  return(X)
}

plot.mobility.margins = function(S,SR){
  S$away.hrs   = S$away.time*24
  S$devices.vt.ht = S$devices.visit / S$devices.home
  S$devices.vt.v0 = S$devices.visit / SR$devices.visit
  S$devices.ht.h0 = S$devices.home  / SR$devices.home
  S$devices.vr.hr = S$devices.vt.v0 / S$devices.ht.h0
  S$devices.vt.h0 = S$devices.visit / SR$devices.home
  # print(summary(S[!idx.ref,]$devices.vr.hr))
  # print(summary(100*S[ idx.ref,]$devices.home / S[ idx.ref,]$devices.total))
  # print(summary(100*S[!idx.ref,]$devices.home / S[!idx.ref,]$devices.total))
  # print(summary(100*(S$pop - S$devices.total) / S$pop))
  # print(summary(S$devices.vt.h0[idx.ref]))
  g = plot.ridge.density(list(
    'Home Reduction'         = data.frame(month=S$month,Ratio=S$devices.ht.h0),
    'Visit Reduction'        = data.frame(month=S$month,Ratio=S$devices.vt.v0),
    'Visit / Reference Home' = data.frame(month=S$month,Ratio=S$devices.vt.h0)
  ),'Ratio',xmax=2); ggsave(figname('devices-panel','mx'),w=15,h=5)
  g = plot.ridge.density(S,x='away.hrs') + labs(x='Hours away from home');
    ggsave(figname('t-away-hours','mx'),w=5,h=5)
  g = plot.ridge.density(S,x='away.ratio') + labs(x='Relative time away from home');
    ggsave(figname('t-away-relative','mx'),w=5,h=5)
}

save.mobility = function(X,S,key='all'){
  # NOTE: this is a bit dangerous since we rely on implicit replication of S
  #       to match X (hence reordering) but it is much faster than merging
  P.g. = aggregate(pop~decile+month,S,sum,na.rm=TRUE,drop=FALSE)$pop
  X = X[order(X$FSA.visited,X$month,X$FSA),]; X$decile = S$decile
  X = X[order(X$FSA,X$month,X$FSA.visited),]; X$decile.visited = S$decile
  X.gg. = aggregate(visit.total~decile+decile.visited+month,X,sum,na.rm=TRUE,drop=FALSE)
  X.gg. = X.gg.[order(X.gg.$decile.visited,X.gg.$month,X.gg.$decile),]
  X.gg.$visit.prop = X.gg.$visit.total / P.g.
  X.gg.$visit.total = NULL
  X.gg. = X.gg.[order(X.gg.$month,X.gg.$decile,X.gg.$decile.visited),]
  write.csv(X.gg.,root.path('data','mix',paste0('mobility_decile_',key,'.csv')),row.names=FALSE)
}

main.mobility = function(do.plot=FALSE,do.csv=TRUE){
  H  = load.fsa.t.away()
  X  = load.fsa.mob()
  S  = merge.mobility.margin(X,H)
  SR = gen.mobility.refs(S,S$month %in% config$t.ref)
  S$away.ratio = S$away.time / SR$away.time
  if (do.csv){
    X = mobility.scale.up(X,S,SR)
    save.mobility(X,S,key='inter')
    X = add.within.fsa.mobility(X,S,SR)
    X = mobility.scale.up(X,S,SR)
    save.mobility(X,S,key='all')
  }
  if (do.plot){
    plot.mobility.margins(S,SR)
  }
}
