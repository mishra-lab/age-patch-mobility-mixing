# supporting functions for loading & cleaning various data sources

aggr.age = function(pop,age=TRUE){
  if (age){
    pop = aggregate(pop~FSA+cut(age,breaks=c(config$age,Inf),right=FALSE),pop,sum)
    colnames(pop)[2] = 'age'
    levels(pop$age) = names(config$age)
    return(pop)
  } else {
    return(aggregate(pop~FSA,pop,sum))
  }
}

aggr.mob.decile = function(B,pop){
  P = aggregate(pop~decile+group,pop,sum)
  p = ave(P$pop,P$group,FUN=function(pop){pop/sum(pop)})
  A = t(tail(outer(c(0,config$group),seq(10),'<=') * outer(c(config$group,Inf),seq(10),'>'),-1))
  B = t(p*A) %*% B %*% A
  colnames(B) = config$group
  rownames(B) = config$group
  return(B)
}

map.decile = function(x){
  x$group = cut(x$decile,breaks=c(config$group,Inf),right=FALSE)
  levels(x$group) = names(config$group)
  return(x)
}

load.fsa.pop = function(age=TRUE,aggr=TRUE){
  dec = read.csv(root.path('data','fsa','fsa_decile.csv'))
  pop = read.csv(root.path('data','fsa','pop_age_fsa.csv'))
  if (aggr){ pop = aggr.age(pop,age=age) }
  pop = merge(pop,map.decile(dec))
  if (age){
    return(pop[order(pop$FSA,pop$age),])
  } else {
    return(pop[order(pop$FSA),])
  }
}

load.group.mob = function(pop){
  mob = read.csv(root.path('data','fsa','mobility_decile.csv'))
  mob = map.decile(mob)
  months = levels(mob$month)
  B.gg.t = lapply(months,function(mo){
    B.m = matrix(mob[mob$month == mo,]$visit.prop,nrow=10,ncol=10)
    return(aggr.mob.decile(B.m,pop))
  })
  names(B.gg.t) = months
  B.gg.t[['REF']] = Reduce('+',B.gg.t[config$t.ref]) / length(config$t.ref)
  B.gg.t = B.gg.t[c('REF',config$t.covid)]
  return(B.gg.t)
}

load.fsa.mob = function(refresh=FALSE){
  rdata = root.path('data','fsa','.rdata','mobility_fsa.rdata')
  if (file.exists(rdata) & !refresh){
    load(rdata)
  } else {
    FSA = levels(load.fsa.pop(age=FALSE)$FSA)
    X = read.csv(root.path('data','fsa','mobility_fsa.csv'))
    colnames(X) = c('FSA.visited','FSA','month','devices.visit','visit.prop','province','devices.home')
    X = X[X$FSA %in% FSA & X$FSA.visited %in% FSA,] # remove external travel
    X$visit.prop = NULL
    X$province = NULL
    X = merge(expand.grid(FSA=FSA,FSA.visited=FSA,month=levels(X$month)),X,all.x=TRUE)
    save(X,file=rdata)
  }
  return(X)
}

load.fsa.smartphones = function(){
  X = read.csv(root.path('data','fsa','smartphones_fsa.csv'))
  return(X)
}

load.contacts = function(c.map=NULL,P.norm=TRUE,sym=TRUE){
  if (is.null(c.map)){ c.map = config$c.map }
  load(root.path('data','fsa','.rdata','Prem2017.rdata'))
  C.AA.y = list()
  for (y in names(c.map)){
    C.AA = Reduce('+',Prem2017$C.AA.y[c.map[[y]]])
    if (P.norm){ C.AA = sweep(C.AA,2,Prem2017$P.a,'/') * mean(Prem2017$P.a) }
    if (sym){    C.AA = symmetric(C.AA) }
    dimnames(C.AA) = list('a'=names(config$age.contact),'a.'=names(config$age.contact))
    C.AA.y[[y]] = C.AA
  }
  names(C.AA.y) = names(config$c.type)
  return(C.AA.y)
}

load.cases = function(){
  X = read.csv(root.path('data','covid','new-cases.csv'))
  X$decile = as.factor(X$decile)
  return(X)
}

clean.raw.pop = function(){
  fsa = read.csv(root.path('data','fsa','fsa.csv'))$FSA
  pop = read.csv(root.path('data','fsa','.raw','pop_age_fsa.csv'))
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
  write.csv(pop,root.path('data','fsa','pop_age_fsa.csv'),row.names=FALSE)
}

clean.Prem2017 = function(){
  # TODO: why does rdata school =/= Prem2017 appendix for a > 30
  rdata = function(name){
    load(root.path('data','fsa','.raw',paste0(name,'.rdata')))
    return(get(name))
  }
  P.raw = rdata('poptotal')
  c.types = list('home','work','school','others')
  names(c.types) = c.types
  Prem2017 = list(
    P.a = as.numeric(P.raw[P.raw$countryname=='Canada',paste0('age',seq(0,75,5))]),
    C.AA.y = lapply(c.types,function(c.type){
      return(rdata(paste0('contact_',c.type))$CAN)
    })
  )
  # plot.mix(Prem2017$C.AA.y,'Ci','a',trans='sqrt'); ggsave('Rplots.pdf',w=14,h=4) # DEBUG
  save(Prem2017,file=root.path('data','fsa','.rdata','Prem2017.rdata')) # TODO: save as CSV?
}