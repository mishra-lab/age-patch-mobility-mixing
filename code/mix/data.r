# supporting functions for loading & cleaning various data sources

aggr.age = function(pop,age=TRUE){
  if (age){
    pop = aggregate(pop~FSA+cut(age,breaks=config$age,right=FALSE),pop,sum)
    colnames(pop)[2] = 'age'
    levels(pop$age) = age.names(config$age)
    return(pop)
  } else {
    return(aggregate(pop~FSA,pop,sum))
  }
}

aggr.mob.decile = function(B,pop){
  P = aggregate(pop~decile+group,pop,sum)
  p = ave(P$pop,P$group,FUN=function(pop){pop/sum(pop)})
  A = t(tail(outer(c(0,config$group),seq(10),'<=') * outer(c(config$group,11),seq(10),'>'),-1))
  B = t(p*A) %*% B %*% A
  dimnames(B)=list(g=config$group,g.=config$group)
  return(B)
}

map.decile = function(x){
  x$group = cut(x$decile,breaks=c(config$group,11),right=FALSE)
  levels(x$group) = names(config$group)
  return(x)
}

load.fsa.pop = function(age=TRUE,aggr=TRUE){
  dec = read.csv(root.path('data','fsa_decile.csv'))
  pop = read.csv(root.path('data','pop_age_fsa.csv'))
  if (aggr){ pop = aggr.age(pop,age=age) }
  pop = merge(pop,map.decile(dec))
  if (age){
    return(pop[order(pop$FSA,pop$age),])
  } else {
    return(pop[order(pop$FSA),])
  }
}

load.decile.pop = function(cat=1){
  pop = read.csv(root.path('data','pop_age_decile.csv'))
  pop = pop[pop$cat==cat,]; pop$cat = NULL
  pop = map.decile(pop)
  return(pop[order(pop$decile,pop$age),])
}

load.group.mob = function(pop,B='Ba',gen=TRUE,rho.src='veraset',
                          months=NULL,do.ref=TRUE,rm.ref=TRUE){
  if (gen){ # generate by approx factors: cond, rel, t
    Xgg   = expand.grid(decile=seq(10),decile.visited=seq(10))
    Bc    = read.csv(root.path('data','mobility_cond.csv'))
    Rg    = read.csv(root.path('data','mobility_rel.csv'))
    rho.t = read.csv(root.path('data','mobility_t.csv'))
    Xgg   = merge(merge(Xgg,Bc),merge(rho.t[rho.t$src==rho.src,],Rg))
    igg = Xgg$decile == Xgg$decile.visited
    Xgg$Ba = Xgg$rho.t * Xgg$R.rho.g *
      ((Xgg$phi.0 * Xgg$R.phi.0.g) * igg + (1 - (Xgg$phi.0 * Xgg$R.phi.0.g)) * Xgg$Bc)
  } else { # observed
    Xgg = read.csv(root.path('data','mobility_obs.csv'))
  }
  Xgg = Xgg[order(Xgg$month,Xgg$decile,Xgg$decile.visited),]
  if (is.null(months)){ months = unique(Xgg$month) }
  B.gg.t = lapply(months,function(month){
    B.gg = t(matrix(Xgg[Xgg$month==month,B],nrow=10,ncol=10))
    dimnames(B.gg) = list(g=seq(10),g.=seq(10))
    return(aggr.mob.decile(B.gg,pop))
  })
  names(B.gg.t) = months
  if (do.ref){ B.gg.t[['REF']] = Reduce('+',B.gg.t[config$t.ref]) / length(config$t.ref) }
  if (rm.ref){ for (month in config$t.ref){ B.gg.t[month] = NULL } }
  return(B.gg.t)
}

load.fsa.mob = function(refresh=FALSE){
  rdata = root.path('data','.rdata','mobility_fsa.rdata')
  if (file.exists(rdata) & !refresh){
    load(rdata)
  } else {
    FSA = levels(load.fsa.pop(age=FALSE)$FSA)
    X = read.csv(root.path('data','.raw','inter_fsa_monthly2.csv'))
    colnames(X) = c('FSA.visited','FSA','month','devices.visit','visit.prop','devices.home')
    X = X[X$FSA %in% FSA & X$FSA.visited %in% FSA,] # remove external travel
    X = X[X$month %in% c(config$t.ref,config$t.covid),]
    X$month = factor(X$month)
    X$visit.prop = NULL
    X = merge(expand.grid(FSA=FSA,FSA.visited=FSA,month=levels(X$month)),X,all.x=TRUE)
    save(X,file=rdata)
  }
  return(X)
}

load.fsa.t.away = function(){
  X = read.csv(root.path('data','.raw','t_away_fsa.csv'))
  X$t.away.inter = X$time.away.inter.mean
  X$t.away.intra = X$time.away.intra.mean
  X$t.away.total = X$t.away.intra + X$t.away.inter
  X[,grepl('time\\.away\\.',names(X))] = NULL
  return(X)
}

load.contacts = function(c.map=NULL,P.norm=TRUE,sym=TRUE){
  if (is.null(c.map)){ c.map = config$c.map }
  load(root.path('data','.rdata','Prem2021.rdata'))
  C.AA.y = list()
  for (y in names(c.map)){
    C.AA = Reduce('+',Prem2021$C.AA.y[c.map[[y]]])
    if (P.norm){ C.AA = sweep(C.AA,2,Prem2021$P.a,'/') * mean(Prem2021$P.a) }
    if (sym){    C.AA = symmetric(C.AA) }
    a.names = age.names(config$age.contact)
    dimnames(C.AA) = list('a'=a.names,'a.'=a.names)
    C.AA.y[[y]] = C.AA
  }
  names(C.AA.y) = names(config$c.type)
  return(C.AA.y)
}

load.cases = function(){
  X = read.csv(root.path('data','covid_cases_decile.csv'))
  X$decile = as.factor(X$decile)
  return(X)
}

load.shape = function(){
  library('sf')
  load(root.path('data','.raw','lfsa000b16a_e.rdata'))
  dec = read.csv(root.path('data','fsa_decile.csv'))
  return(merge(shp,dec,all.y=TRUE))
}

clean.raw.pop = function(){
  fsa = read.csv(root.path('data','fsa.csv'))$FSA
  pop = read.csv(root.path('data','pop_age_fsa.csv')) # TODO: this was .raw before, why?
  pop$age = as.character(pop$Age_cat)
  pop$Age = NULL
  pop$Age_cat = NULL
  pop$CENSUS_YEAR = NULL
  pop$Female = NULL
  pop$Male = NULL
  pop$age[pop$age == 'Under 1 year'] = '0'
  pop$age[pop$age == '100 years and over'] = '100'
  pop = pop[grepl('^\\d+$',pop$age),]
  pop$age = as.numeric(pop$age)
  pop = col.rename(pop,'Total','pop')
  pop = col.rename(pop,'GEO_NAME','FSA')
  pop = pop[pop$FSA %in% fsa,]
  write.csv(pop,root.path('data','pop_age_fsa.csv'),row.names=FALSE)
}

clean.Prem2021 = function(){
  rdata = function(name){
    load(root.path('data','.raw',paste0(name,'.rdata')))
    return(get(name))
  }
  P.raw = rdata('poptotal')
  c.types = list('home','work','school','others')
  names(c.types) = c.types
  Prem2021 = list(
    P.a = as.numeric(P.raw[P.raw$countryname=='Canada',paste0('age',seq(0,75,5))]),
    C.AA.y = lapply(c.types,function(c.type){
      return(rdata(paste0('contact_',c.type))$CAN)
    })
  )
  # plot.mix(Prem2021$C.AA.y,'Ci','a',trans='sqrt'); ggsave('Rplots.pdf',w=14,h=4) # DEBUG
  save(Prem2021,file=root.path('data','.rdata','Prem2021.rdata'))
}

clean.shape = function(which='FSA'){
  fname = function(ext){ root.path('data','.raw',paste0('lfsa000b16a_e.',ext)) }
  shp = sf::st_read(dsn=fname('shp'),quiet=TRUE)
  shp = subset(shp,PRNAME=='Ontario')
  shp$FSA = shp$CFSAUID
  shp$CFSAUID = NULL
  shp$PRNAME = NULL
  shp$PRUID = NULL
  save(shp,file=fname('rdata'))
}