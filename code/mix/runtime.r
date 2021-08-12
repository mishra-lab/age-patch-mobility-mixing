# Objective: Same as mixing.r except pre-compute as much as possible
#            and save/load things from (fsa.mix.data(mode)) .rdata file
# 1. Run "gen.mix.runtime.data" once to generate the .rdata file
# 2. Run "gen.Ci.gaga.y" during model fitting with input values sampled as needed


gen.mix.runtime.data = function(){
  # Pre-compute as much mixing data as possible and save in .rdata file
  source('mixing.r')
  config = set.config(mode='10x10',n.y='2')
  pop    = load.fsa.pop()
  P.ga   = pop.to.Pga(pop)
  C.AA.y = load.contacts()
  C.aa.y = CAAy.to.Caay(C.AA.y)
  a.size = bin.size(config$age,norm=TRUE)
  B.gg.t = load.group.mob(pop)
  save(P.ga,C.aa.y,a.size,B.gg.t,config,file='.fsa-mix.rdata')
}

a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}

gen.mix.runtime = function(t='REF',RC.global=NULL,RC.decile=NULL){
  load('.fsa-mix.rdata')
  if (missing(RC.global)){ RC.global = config$RC.global }
  if (missing(RC.decile)){ RC.decile = config$RC.decile }
  N = config$N; h.y = config$h.y; # convenience
  B.gg = B.gg.t[[t]] # mobility during t
  B.hh = diag(N$g) # dummy "mobillity" for household contacts
  X.gaga.y = list()
  for (y in seq(N$y)){ # for each contact type
    X.gaga = array(0,c(N$g,N$a,N$g,N$a),dimnames=config$X.names) # initialize X for this y
    C.gaga = aperm(replicate(N$g,replicate(N$g,C.aa.y[[y]])),c(3,1,4,2)) # pad C.aa to C.gaga
    RC.g   = RC.global[y] * RC.decile/mean(RC.decile)
    for (g in seq(N$g)){ # for each group
      # traveller pool
      Pt  = RC.g * P.ga  * (1-h.y[y]) * B.gg[,g] # pop of each "g" mixing in this pool
      pt. = sweep(Pt,2,a.size,'/')/sum(Pt) # normalize overall & age as C.gaga already weighted by age
      pt.[is.nan(pt.)] = 0
      X = outer(Pt,pt.) * C.gaga # proportionate mixing * age mixing
      X.gaga = X.gaga + X # add contribution of this pool to the total mixing
      # home pool
      Ph  = RC.g * P.ga[g,]  * h.y[y] # pop (age groups) of this "g" at home
      ph. = Ph/a.size/sum(Ph) # normalize overall & age as C.aa.y already weighted by age
      ph.[is.nan(ph.)] = 0
      X = outer(Ph,ph.) * C.aa.y[[y]] # proportionate mixing * age mixing
      X.gaga[g,,g,] = X.gaga[g,,g,] + X # add contribution of this pool to the total mixing
    }
    X.gaga = X.gaga/2 + aperm(X.gaga,c(3,4,1,2))/2 # ensure contacts balance
    for (g in seq(N$g)){
      for (a in seq(N$a)){
        X.gaga[g,a,,] = X.gaga[g,a,,] / P.ga[g,a]
      }
    }
    X.gaga.y[[y]] = X.gaga # copy result into output
  }
  names(X.gaga.y) = names(config$c.type)
  return(X.gaga.y) # "CX"
}

gen.mix.original = function(){
  source('mixing.r')
  return(gen.save.mixing(do.save=FALSE)$REF)
}

test.runtime = function(){
  gen.mix.runtime.data()
  Ci.runtime = gen.mix.runtime()
  Ci.original = gen.mix.original()
  for (y in seq(config$N$y)){ # only print if we fail
    if (isTRUE(all.equal(Ci.original[[y]],Ci.runtime[[y]]))){
      print('Ci.original == Ci.runtime [OK]')
    } else {
      print(paste('Ci.original =/= Ci.runtime: y =',y))
    }
  }
}