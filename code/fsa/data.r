source('config.r')

load.fsa.pop = function(age=TRUE){
  pop = read.csv(root.path('data','fsa','age_pop.csv'))
  if (age){
    return(pop)
  } else {
    FSA = list(FSA=sort(pop$FSA))
    return(merge(
      aggregate(list(pop=pop$pop_agefsa),FSA,sum,drop=FALSE),
      aggregate(list(decile=pop$decile),FSA,mean,drop=FALSE)
    ))
  }
}

load.decile.mob = function(){
  mob = read.csv(root.path('data','fsa','mobility_decile.csv'))
  months = levels(mob$month)
  B.gg = lapply(months,function(mo){
    B.m = matrix(mob[mob$month == mo,]$visit.prop,nrow=10,ncol=10)
    colnames(B.m) = info$decile
    rownames(B.m) = info$decile
    return(B.m)
  })
  names(B.gg) = months
  return(B.gg)
}

load.fsa.mob = function(refresh=FALSE){
  rdata = root.path('data','fsa','.rdata','mobility.fsa.rdata')
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
  X = read.csv(root.path('data','fsa','smartphones.csv'))
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
