source('config.r')

aggr.age = function(pop,age=TRUE){
  if (age){
    pop = aggregate(pop~FSA+cut(age,breaks=c(info$age,Inf),right=FALSE),pop,sum)
    colnames(pop)[2] = 'age'
    levels(pop$age) = names(info$age)
    return(pop)
  } else {
    return(aggregate(pop~FSA,pop,sum))
  }
}

aggr.mob.decile = function(B,pop){
  P = aggregate(pop~decile+group,pop,sum)
  p = ave(P$pop,P$group,FUN=function(pop){pop/sum(pop)})
  M = t(tail(outer(c(0,info$group),seq(10),'<=') * outer(c(info$group,Inf),seq(10),'>'),-1))
  B = t(p*M) %*% B %*% M
  colnames(B) = info$group
  rownames(B) = info$group
  return(B)
}

map.decile = function(x){
  x$group = cut(x$decile,breaks=c(info$group,Inf),right=FALSE)
  levels(x$group) = names(info$group)
  return(x)
}

load.fsa.pop = function(age=TRUE){
  dec = read.csv(root.path('data','fsa','fsa_decile.csv'))
  pop = read.csv(root.path('data','fsa','pop_age_fsa.csv'))
  pop = merge(aggr.age(pop,age=age),map.decile(dec))
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
  B.gg = lapply(months,function(mo){
    B.m = matrix(mob[mob$month == mo,]$visit.prop,nrow=10,ncol=10)
    return(aggr.mob.decile(B.m,pop))
  })
  names(B.gg) = months
  return(B.gg)
}

load.fsa.mob = function(refresh=FALSE){
  rdata = root.path('data','fsa','.rdata','mobility_fsa.rdata')
  if (file.exists(rdata) & !refresh){
    load(rdata)
  } else {
    X = read.csv(root.path('data','fsa','mobility_fsa.csv'))
    X = X[X$visited_fsa %in% unique(X$home_fsa),] # remove external travel
    colnames(X) = c('FSA.visited','FSA','month','devices.visit','visit.prop','devices.home')
    X$visit.prop = NULL
    FSA   = sort(unique(X$FSA))
    month = sort(unique(X$month))
    X = merge(expand.grid(FSA=FSA,FSA.visited=FSA,month=month),X,all=TRUE)
    save(X,file=rdata)
  }
  return(X)
}

load.fsa.smartphones = function(){
  X = read.csv(root.path('data','fsa','smartphones_fsa.csv'))
  return(X)
}

load.contacts = function(){
  C.y = list()
  map = list(Home='home',Other=c('work','school','other_locations'))
  ctx = read.csv(root.path('data','fsa','contacts_age.csv'))
  dnames = list(a=names(age.contact),a.=names(age.contact))
  for (y in names(map)){
    C.y[[y]] = matrix(
      aggregate(contacts~a+a.,ctx[ctx$type %in% map[[y]],],sum)$contacts,
      nrow=length(age.contact),dimnames=dnames)
  }
  # g = plot.mix(C.y,aggr=FALSE); ggsave(figname('polymod-canada-2'),width=8,height=4) # DEBUG
  return(C.y)
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

clean.canada.contacts = function(){
  library('readxl')
  library('reshape2')
  C.y = list()
  for (y in c('home','work','school','other_locations')){
    f = root.path('data','fsa','.raw','contacts-152',paste0('MUestimates_',y,'_1.xlsx'))
    C.y. = as.matrix(read_excel(f,sheet='Canada'))
    rownames(C.y.) = names(age.contact)
    colnames(C.y.) = names(age.contact)
    # C.y[[y]] = C.y. # DEBUG
    C.y[[y]] = cbind('type'=y,melt(C.y.,value.name='contacts',varnames=c('a','a.')))
  }
  # C.y. = list(Home=C.y[['home']],Other=C.y[['work']]+C.y[['school']]+C.y[['other_locations']]) # DEBUG
  # plot.mix(C.y,aggr=FALSE,clim=c(0,7)); ggsave(figname('polymod-canada'),width=12,height=4) # DEBUG
  C.y = do.call(rbind,C.y)
  write.csv(C.y,root.path('data','fsa','contacts_age.csv'),row.names=FALSE)
}