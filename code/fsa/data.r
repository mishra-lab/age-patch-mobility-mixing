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

load.polymod = function(){
  pmc = list()
  f = 1
  map = list(Home=c('home'),Other=c('work','school','transport','leisure','otherplace'))
  for (c.type in info$c.type){
    pmc[[c.type]] = matrix(0,nrow=length(info$age),ncol=length(info$age))
    rownames(pmc[[c.type]]) = names(info$age)
    colnames(pmc[[c.type]]) = names(info$age)
    for (type in map[[c.type]]){
      names(f) = paste0('cnt_',type)
      pmc[[c.type]] = pmc[[c.type]] +
        suppressMessages(contact_matrix(polymod,filter=f,age.limits=info$age))$matrix
    }
  }
  return(pmc)
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
