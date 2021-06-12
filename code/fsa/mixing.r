source('utils.r')
source('data.r')
source('plot.r')

pop.to.Pga = function(pop){
  # pop: dataframe with columns (at least): group, age, pop(ulation)
  # P.ga: matrix of population by group & age
  P = aggregate(pop~group+age,pop,sum)
  P.ga = matrix(P[order(P$age,P$group),]$pop,nrow=config$N$g,
    dimnames=list(g=names(config$group),a=names(config$age)))
  return(P.ga)
}

Caay.to.Cay.resample = function(C.aa.y){
  # C.aa.y: list by contact type of matrices: contacts by age & other age
  #         but not necessarily the right age bins
  # C.a.y:  list by contact type of vectors: contacts by age, with the right bins
  midpoint = function(a){ a+c(diff(a)/2,0) }
  C.a.y = lapply(C.aa.y,function(C.aa){ # for each contact type
    c.a.i = rowSums(C.aa) # total contacts by age group (sum across "other" ages)
    c.a.o = approx( # linearly interpolate c.a.i (original age bins) to get c.a.o (new age bins)
      y  = c.a.i,
      x  = midpoint(config$age.contact),
      xo = midpoint(config$age),
      method = 'linear',
      yleft  = c.a.i[1],
      yright = c.a.i[length(c.a.i)])$y
    names(c.a.o) = names(config$age)
    return(c.a.o)
  })
  return(C.a.y)
}

Cay.to.Cgay = function(C.a.y){
  # C.a.y:  list by contact type of vectors: contacts by age
  # C.ga.y: list by contact type of matrices: contacts by group & age
  S = map.decile(data.frame(decile=seq(10)))
  C.ga.y = lapply(names(C.a.y),function(y){ # for each contact type
    S$scale = config$f.c.mean[[y]]+seq(+config$f.c.slope[[y]],-config$f.c.slope[[y]],l=10) # C scaling
    df = transform(merge(data.frame(C=C.a.y[[y]],age=config$age),S),CS=C*scale) # CS = scaled contacts
    CS = aggregate(CS~group+age,df,mean)$CS # aggregate CS by group (easier to do after scaling)
    return(matrix(CS,nrow=config$N$g,dimnames=config$X.names[1:2])) # reshape as matrix by group & age
  })
  names(C.ga.y) = names(C.a.y)
  return(C.ga.y)
}

gen.mix.total = function(C.ga.y,P.ga,B.gg.t,t){
  # estimates the total number of contacts between all combinations of age/geo groups
  # C.ga.y: list by contact type of matrices: contacts by group & age
  # P.ga:   matrix of population by group & age
  # B.gg.t: list by time of matrices: mobility by group & other group
  # t:      month
  # X.gaga.y: list by type of matrices: total # contacts by age, group, other age, other group (not per person)
  B.gg.t[['REF']] = Reduce('+',B.gg.t[config$t.ref]) / length(config$t.ref) # mobility in REF
  B.g.0 = rowSums(B.gg.t[['REF']]) # total probability of being mobile during REF
  B.g.t = rowSums(B.gg.t[[t]])     # total probability of begin mobile this t
  B.gg = B.gg.t[[t]] + (B.g.t / B.g.0) * diag(1-B.g.0) # adding "mobility" within same g
  Q.ga.y = lapply(config$c.type,function(y){ P.ga * C.ga.y[[y]] }) # total contacts offered by group & age, by type
  N = config$N; eps.y = config$eps.y; h.y = config$h.y; # convenience
  X.gaga.y = list() # initialize output
  for (y in seq(N$y)){ # for each contact type
    X.gaga = array(0,c(N$g,N$a,N$g,N$a),dimnames=config$X.names) # initialize output for thie type
    # non-home contacts
    for (. in seq(N$g)){ # for each mixing pool (group)
      Q.ga    = Q.ga.y[[y]] * B.gg[,.] # contacts made available to the pool from all groups by group & age
      X.gaga. = outer(Q.ga,Q.ga) / sum(Q.ga) # proportionate distribution of contacts
      X.gag   = a.sum(X.gaga.,4) # sum across "other" age groups
      X.gaga. = X.gaga. * (1-eps.y[[y]]) # proportionate part of age-mixing for this pool
      for (a in seq(N$a)){ # for each age group
        X.gaga.[,a,,a] = X.gaga.[,a,,a] + eps.y[[y]] * X.gag[,a,] # assortative part of age-mixing for this pool
      }
      X.gaga = X.gaga + (1-h.y[[y]]) * X.gaga. # add contribution of this away-pool to the total mixing
    }
    # home contacts
    for (g in seq(N$g)){ # for each group
      Q.a   = Q.ga.y[[y]][g,] # contacts made available to this pool (no mobility/decile involved)
      X.aa. = outer(Q.a,Q.a) / sum(Q.a) # proportionate distribution of contacts
      X.a   = a.sum(X.aa.,2) # sum across "other" age groups
      X.aa. = X.aa. * (1-eps.y[[y]]) # proportionate part of age-mixing for this pool
      for (a in seq(N$a)){ # for each age group
        X.aa.[a,a] = X.aa.[a,a] + eps.y[[y]] * X.a[a] # assortative part of age-mixing for this pool
      }
      X.gaga[g,,g,] = X.gaga[g,,g,] + (h.y[[y]]) * X.aa. # add contribution of this home-pool to the total mixing
    }
    X.gaga.y[[y]] = X.gaga # copy result into output 
  }
  names(X.gaga.y) = names(config$c.type)
  # for (y in seq(N$y)){ print(a.sum(X.gaga.y[[y]],c(3,4))/Q.ga.y[[y]]) } # DEBUG: should all = 1
  return(X.gaga.y) # "CX"
}

CX.norm = function(CX.gaga.y,P.ga,C.ga.y=NULL){
  # convert CX (total contacts) to:
  # - Ci: per person, divide by P.ga (if C.ga.y is NULL)
  # - Cp probability, divide by P.ga & C.ga (if C.ga.y is given)
  C.gaga.y = lapply(CX.gaga.y,function(X){ X*0 }) # initialize output
  for (y in seq(config$N$y)){
    for (g in seq(config$N$g)){
      for (a in seq(config$N$a)){
        c.ga.y = ifelse(is.null(C.ga.y),1,C.ga.y[[y]][g,a])
        C.gaga.y[[y]][g,a,,] = CX.gaga.y[[y]][g,a,,] / P.ga[g,a] / c.ga.y
      }
    }
  }
  names(C.gaga.y) = names(CX.gaga.y)
  return(C.gaga.y)
}

aggr.mix = function(C.gaga,what,vs,P.ga=NULL,aggr=TRUE){
  # aggregate C.gaga (4D) to produce:
  # - what='a': age vs other age (aggregate across group & other group)
  # - what='g': group vs other group (aggregate across age & other age)
  # - what='i': group vs age (aggregate across other group & other age)
  # - what='o': othe group vs other age (not implemented, not sure the interpretation)
  # aggregation operation depends on if C.gaga is "CX", "Ci", "Cp"
  if (!aggr){ return(C.gaga) }
  d = switch(vs,'a'=c(1,3),'g'=c(2,4),'i'=c(3,4)) # 'o'=c(1,2) TODO: check interpretation of 'o'
  if (what=='CX'){ return( a.sum(C.gaga,d) / 1e6 )  } # sum absolute contacts & divide by 1 million
  if (what=='Ci'){ # sum across "other" dimensions, weighted average across "self" dimensions
    P.gaga = array(P.ga,dim(C.gaga))
    return( a.mean( a.sum(C.gaga,d[d>=3]), d[d<=2], a.sum(P.gaga,d[d>=3])) )
  }
  if (what=='Cp'){ stop('aggr.mix for Cp not yet implemented') }
}

prop.self = function(X){ # what proportion of X is along the diagonal?
  return(median(diag(X) / rowSums(X)))
}

plot.mixing = function(C.gaga.y,what,t='REF',...){
  # print(sapply(C.gaga.y,function(C){ prop.self( aggr.mix(C,what,'g',...) ) }))
  plot.mix(C.gaga.y,what,'g',...);            ggsave(figname(paste0(what,'ggy'),  'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'g',...,xfun=offd);  ggsave(figname(paste0(what,'ggy-o'),'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'a',...);            ggsave(figname(paste0(what,'aay'),  'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'a',...,xfun=offd);  ggsave(figname(paste0(what,'aay-o'),'mixing',MODE,t),width=8,height=4)
}

mixing.fname = function(what,t,sub='.tmp/mix'){
  # standard filename for saving data
  path = root.path('data','fsa',sub)
  suppressWarnings({ dir.create(path,recursive=TRUE) })
  return(file.path(path,paste0(what,'_',MODE,'_',t,'.csv')))
}

save.mixing = function(X.gaga.y,what,t){
  # melt the contents of X.gaga.y to csv with colnames shown below
  mix. = do.call(expand.grid,dimnames(X.gaga.y[[1]]))
  mix = do.call(rbind,lapply(names(X.gaga.y),function(y){
    cbind(mix.,y=config$c.type[[y]],X=as.vector(X.gaga.y[[y]]))
  }))
  v.name = switch(what,CX='n.contacts',Ci='n.contact.pp',Cp='prob.contact')
  colnames(mix) = c('index.decile','index.age','other.decile','other.age','contact.type',v.name)
  write.csv(mix,mixing.fname(what,t),row.names=FALSE)
}

merge.save.mixing = function(what){
  # main.mixing only does one t; this function merges the csv files & re-saves as one
  mix = do.call(rbind,lapply(c('REF',config$t.covid),function(t){
    mix. = cbind(month=t,read.csv(mixing.fname(what,t)))
  }))
  write.csv(mix,mixing.fname(what,'all',''),row.names=FALSE)
}

main.mixing = function(t='REF'){
  pop    = load.fsa.pop()
  C.aa.y = load.contacts()
  B.gg.t = load.group.mob(pop)
  C.ga.y = Cay.to.Cgay(Caay.to.Cay.resample(C.aa.y))
  P.ga   = pop.to.Pga(pop)
  CX.gaga.y = gen.mix.total(C.ga.y,P.ga,B.gg.t,t) # total absolute contacts
  Ci.gaga.y = CX.norm(CX.gaga.y,P.ga)             # contacts per person
  # Cp.gaga.y = CX.norm(CX.gaga.y,P.ga,C.ga.y)      # contact probability
  plot.mixing(CX.gaga.y,'CX',t)
  plot.mixing(Ci.gaga.y,'Ci',t,P.ga=P.ga)
  save.mixing(Ci.gaga.y,'Ci',t)
}
