source('data.r')
source('plot.r')

check.align = function(X1,X2,...){ # for DEBUG only
  args = list(...)
  for (arg in args){
    print(mean(X1[[arg]]==X2[[arg]]))
  }
}

CoV.by = function(X,name,by){
  X.by = lapply(by,function(i){X[[i]]})
  print(summary(sapply(split(X,X.by),function(Xi){
    sd(Xi[[name]]) / mean(Xi[[name]])
  })))
}

plot.mobility = function(Tg,Tn,Xgg,Xnn){
  Tgc = Tg[Tg$is.covid,]
  # CoV stats
  CoV.by(Xgg,'Bc',c('decile','decile.visited'))
  CoV.by(Tgc,'R.rho.gt','decile')
  CoV.by(Tg,'R.phi.gt','decile')
  CoV.by(Tg,'R.phi.0.gt','decile')
  # Xgg
  plot.B.approx(Xgg); ggsave(figname('B-vs-Ba'),w=10,h=5)
  plot.cond.mob(Xgg,'line'); ggsave(figname('Bc-t'),w=10,h=5.5)
  plot.cond.mob(Xgg,'box');  ggsave(figname('Bc-box'),w=10,h=5)
  # Tg: p.mob (rho)
  plot.p.mob(Tg,y='rho.t',style='mean');        ggsave(figname('rho-t'),w=3.5,h=3)
  plot.p.mob(Tg,y='rho.gt/rho.t',style='line'); ggsave(figname('R-rho-gt'),w=4,h=3)
  plot.p.mob(Tgc,y='rho.gt/rho.t',style='box'); ggsave(figname('R-rho-box'),w=4,h=3)
  # Tg: p.intra (phi)
  plot.p.mob(Tg,y='phi.t',style='mean');        ggsave(figname('phi-t'),w=3.5,h=3)
  plot.p.mob(Tg,y='phi.gt/phi.t',style='line'); ggsave(figname('R-phi-gt'),w=4,h=3)
  plot.p.mob(Tg,y='phi.gt/phi.t',style='box');  ggsave(figname('R-phi-box'),w=4,h=3)
  # Tn
  plot.ridge.density(Tn,'rho.gt',fill='decile',xmax=2,legend=TRUE) +
    xlab('Relative time away from home')
    ggsave(figname('rho-gnt'),w=5.5,h=5)
  plot.ridge.density(Tn,'phi.gt',fill='decile',xmax=1,legend=TRUE) +
    xlab('Proportion of away time within home FSA')
    ggsave(figname('phi-gnt'),w=5.5,h=5)
  # Xnn
  Xn = gen.device.ratios(Xnn)
  plot.ridge.density(Xn,'R.home',xmax=2) +
    xlab('Reduction in observed devices') + theme(legend.position='top')
    ggsave(figname('RHt'),w=5,h=5)
  plot.ridge.density(Xn,'R.visit',xmax=2) +
    xlab('Reduction in device visits') + theme(legend.position='top')
    ggsave(figname('RVt'),w=5,h=5)
}

gen.device.ratios = function(Xnn){
  Xn = merge(
    aggregate(devices.home~FSA+month,Xnn,mean),
    aggregate(devices.visit~FSA+month,Xnn,sum))
  Xn = Xn[order(Xn$month,Xn$FSA),]
  Xn.ref = aggregate(cbind(devices.home,devices.visit)~FSA,
    Xn[Xn$month %in% config$t.ref,],mean)
  check.align(Xn,Xn.ref,'FSA') # DEBUG == 1
  Xn$R.home  = aggregate(devices.home~FSA+month,Xn,mean)$devices.home / Xn.ref$devices.home
  Xn$R.visit = aggregate(devices.visit~FSA+month,Xn,sum)$devices.visit / Xn.ref$devices.visit
  return(Xn)
}

gen.p.mob = function(T,by,do.extra=TRUE,rm.src=TRUE){
  T$is.covid = T$month %in% config$t.covid
  Ti = aggregate(formula(paste('cbind(t.away.total,t.away.intra) ~',by)),
    T[!T$is.covid,],mean,na.rm=FALSE,drop=TRUE)
  T = T[order(T[['month']],T[[by]]),]
  # check.align(T,Ti,by) # DEBUG == 1
  T$rho.gt   = T$t.away.total / Ti$t.away.total
  T$phi.gt   = T$t.away.intra / T$t.away.total
  if (do.extra){
    T = T[order(T[[by]],T[['month']]),]
    T$rho.t   = aggregate(rho.gt~month,T,mean)$rho.gt
    T$phi.t   = aggregate(phi.gt~month,T,mean)$phi.gt
    T$phi.0   = mean(T$phi.t)
    T = T[order(T[['month']],T[[by]]),]
    T$R.rho.gt = T$rho.gt / T$rho.t
    T$R.phi.gt = T$phi.gt / T$phi.t
    T$R.phi.0.gt = T$phi.gt / T$phi.0
    T$R.rho.g = aggregate(R.rho.gt~decile,T,mean)[[2]]
    T$R.phi.g = aggregate(R.phi.gt~decile,T,mean)[[2]]
    T$R.phi.0.g = aggregate(R.phi.0.gt~decile,T,mean)[[2]]
  }
  if (rm.src){
    T$t.away.total = NULL
    T$t.away.intra = NULL
  }
  return(T)
}

merge.Xnn.pop = function(Xnn,pop){
  Xn = aggregate(devices.visit~FSA+month,Xnn,sum,drop=FALSE)
  Xn = merge(Xn,pop,all.x=TRUE)
  Xn  = Xn[order(Xn$month,Xn$FSA),]
  Xnn = Xnn[order(Xnn$FSA.visited,Xnn$month,Xnn$FSA),];
  Xnn$decile = Xn$decile
  Xnn = Xnn[order(Xnn$FSA,Xnn$month,Xnn$FSA.visited),];
  Xnn$decile.visited = Xn$decile
  return(Xnn)
}

gen.Bc = function(X,by){
  by.v = paste0(by,'.visited')
  X.   = aggregate(formula(paste('devices.visit ~',by,'+ month')),X,sum,drop=FALSE)
  X    = X[order(X[[by.v]],X[['month']],X[[by]]),]
  # check.align(X,X.,'month',by) # DEBUG == 1
  X$Bc = X$devices.visit / X.$devices.visit
  # print(aggregate(formula(paste('Bc ~',by,'+ month')),X,sum)) # DEBUG == 1
  return(X)
}

main.mobility = function(do.plot=TRUE,do.csv=TRUE){
  pop = load.fsa.pop(age=FALSE)
  pop$decile = as.factor(pop$decile)
  # time away data -> relative mobility vs REF (p.mob) & % mobility intra (p.intra)
  Tn = merge(load.fsa.t.away(),pop)
  Tg = aggregate(cbind(t.away.total,t.away.intra)~decile+month,Tn,mean)
  Tg = gen.p.mob(Tg,'decile')
  # inter FSA mobility data -> Bc
  Xnn = load.fsa.mob()
  Xnn = merge.Xnn.pop(Xnn,pop)
  Xgg = aggregate(devices.visit~decile+decile.visited+month,Xnn,sum,na.rm=TRUE)
  Xgg = gen.Bc(Xgg,'decile')
  Xgg = merge(Xgg,Tg)
  igg = Xgg$decile == Xgg$decile.visited
  Xgg$B  = Xgg$rho.gt * (Xgg$phi.gt * igg + (1 - Xgg$phi.gt) * Xgg$Bc)
  Xgg$Ba = Xgg$rho.t * Xgg$R.rho.g *
    ((Xgg$phi.0 * Xgg$R.phi.0.g) * igg + (1 - (Xgg$phi.0 * Xgg$R.phi.0.g)) * Xgg$Bc)
  Xgg = Xgg[order(Xgg$month,Xgg$decile,Xgg$decile.visited),]
  if (do.plot){
    Tn  = gen.p.mob(Tn,'FSA',do.extra=FALSE,rm.src=FALSE)
    # Xnn = gen.Bc(Xnn,'FSA') # unused & expensive
    plot.mobility(Tg,Tn,Xgg,Xnn)
  }
  if (do.csv){
    Xgg = Xgg[,c('decile','month','decile.visited','rho.t','phi.0','Bc','B','Ba')]
    write.csv(Xgg,root.path('data','mix','mobility_decile.csv'),row.names=FALSE)
  }
}