fsa.mix.data = function(mode='10x10'){
  return(paste0('.fsa-mix-',mode,'.rdata'))
}

gen.mix.runtime.data = function(){
  # Pre-compute as much mixing data as possible and save in fsa.mix.data() file
  source('data.r')
  source('mixing.r')
  pop    = load.fsa.pop()
  C.aa.y = load.contacts()
  B.gg.t = load.group.mob(pop)
  P.ga   = pop.to.Pga(pop)
  C.a.y  = Caay.to.Cay.resample(C.aa.y)
  C.df.y = lapply(seq(config$N$y),function(y){merge( # 1st 1/2 of Cay.to.Cgay
    map.decile(data.frame(decile=seq(10))),
    data.frame(C=C.a.y[[y]],age=config$age)
  )})
  B.gg.t[['REF']] = Reduce('+',B.gg.t[config$t.ref]) / length(config$t.ref)
  B.gg.t = B.gg.t[c('REF',config$t.covid)]
  save(B.gg.t,P.ga,C.df.y,config,file=fsa.mix.data(MODE))
}

a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}

gen.Ci.gg.y = function(t='REF', c.mean=c(1,1), c.slope=c(.15,.30), mode='10x10'){
  load(fsa.mix.data(mode))
  N = config$N; eps.y = config$eps.y; h.y = config$h.y; # convenience
  C.ga.y = lapply(seq(N$y),function(y){ # 2nd 1/2 of Cay.to.Cgay
    C.df.y[[y]]$CS = C.df.y[[y]]$C * (c.mean[y]+seq(+c.slope[y],-c.slope[y],l=10)) # C scaling
    CS = aggregate(CS~group+age,C.df.y[[y]],mean)$CS # aggregate CS by group
    return(matrix(CS,nrow=config$N$g,dimnames=config$X.names[1:2])) # reshape as matrix by group & age
  })
  B.g.0 = rowSums(B.gg.t[['REF']]) # total probability of being mobile during REF
  B.g.t = rowSums(B.gg.t[[t]])     # total probability of begin mobile this t
  B.gg = B.gg.t[[t]] + (B.g.t / B.g.0) * diag(1-B.g.0) # adding "mobility" within same g
  Q.ga.y = lapply(seq(N$y),function(y){ P.ga * C.ga.y[[y]] }) # total contacts offered by group & age, by type
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
    # normalize to per-person contacts
    for (g in seq(N$g)){
      for (a in seq(N$a)){
        X.gaga[g,a,,] = X.gaga[g,a,,] / P.ga[g,a]
      }
    }
    X.gaga.y[[y]] = X.gaga # normalize & copy result into output
  }
  names(X.gaga.y) = names(config$c.type)
  return(X.gaga.y) # "Ci"
}