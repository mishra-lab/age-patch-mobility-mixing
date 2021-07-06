source('utils.r')
source('data.r')
source('plot.r')

# Objective: produce age & group (decile) mixing matrix,
#            stratified by contact type, for a specific time period (month)
# Notation: C: contact numbers. C & Ci: per-person; CX: total number in model
#           P: population size
#           B: mobility matrix
# Most things are matrices, where the subscripts are denoted like:
#   *.g: stratified by group = decile if mode = 10x10, or deciles 1-2/3-10 if mode = 2x2
#   *.a: stratified by age
#   *.gg/.aa: stratified by self & other group / age
# A few things are lists of matrices, one for each contact type (*.**.y)
# Everything here runs for one time period (t) only.

pop.to.Pga = function(pop){
  # pop: dataframe with columns (at least): group, age, pop(ulation)
  # P.ga: [value] matrix of population by group & age
  P = aggregate(pop~group+age,pop,sum)
  P.ga = matrix(P[order(P$age,P$group),]$pop,nrow=config$N$g,
    dimnames=list(g=names(config$group),a=names(config$age)))
  return(P.ga)
}

resample.contacts = function(Ci,xi,xo,x.max){
  # Ci:    contact matrix stratified by xi & xi (square)
  # xi:    list of cut points in Ci (lower cut points of each strata)
  # xo:    list of cut points in Co (lower cut points of each strata)
  # x.max: max cut point for both
  # Co:    [value] contact matrix re-stratified by xo * xo (square)
  # We cannot simply interp2d, since the number of contacts per person would change;
  # that's why we interpolate to m. * m. first, then aggregate using M
  # and normalize by the strata sizes (diff);
  # the number of contacts per person might change a little.
  m. = midpoint(seq(0,x.max-1,1))
  C. = interp2d(Ci/bin.size(xi,x.max),midpoint(xi),m.)
  M = tail(outer(c(0,xo),m.,'<') * outer(c(xo,Inf),m.,'>'),-1)
  Co = (M %*% C. %*% t(M)) / bin.size(xo,x.max)
  return(Co)
}

CAAy.to.Caay = function(C.AA.y){
  # C.AA.y: list of contact matrices in original age stratification
  # C.aa.y: [value] list of contact matrices in desired age stratification
  return(lapply(C.AA.y,resample.contacts,
    xi=config$age.contact,
    xo=config$age,
    x.max=config$age.max))
}

gen.mix.main = function(C.aa.y,P.ga,B.gg.t,t){
  # estimates the total number of contacts between all combinations of age/geo groups
  # C.ga.y: list by contact type of matrices: contacts by group & age
  # P.ga:   matrix of population by group & age
  # B.gg.t: list by time of matrices: mobility by group & other group
  # t:      month
  # X.gaga.y: [value] list by type of arrays: total number of contacts
  #           by age, group, other age, other group (not per person)
  N = config$N; h.y = config$h.y; # convenience
  B.g.0 = rowSums(B.gg.t[['REF']]) # total mobility during REF
  B.g.t = rowSums(B.gg.t[[t]])     # total mobility during t
  B.gg = B.gg.t[[t]] + (B.g.t / B.g.0) * diag(1-B.g.0) # add mobility in same "g"
  B.hh = diag(N$g) # dummy "mobillity" for household contacts
  a.size = bin.size(config$age,config$age.max,norm=TRUE) # expected age distribution
  X.gaga.y = list() # initialize output
  for (y in seq(N$y)){ # for each contact type
    X.gaga = array(0,c(N$g,N$a,N$g,N$a),dimnames=config$X.names) # initialize X for this y
    C.gaga = aperm(replicate(N$g,replicate(N$g,C.aa.y[[y]])),c(3,1,4,2)) # pad C.aa to C.gaga
    for (. in seq(N$g)){ # for each mixing pool (group)
      P. = P.ga * ((1 - h.y[y]) * B.gg[,.] + (h.y[y] * B.hh[,.])) # pop of each "g" mixing here
      PR = sweep(P.,2,a.size,'/') / sum(P.) # normalize by expected age distrib & overall pop
      X = outer(P.,PR) * C.gaga # proportionate mixing * age mixing
      X.gaga = X.gaga + X # add contribution of this pool to the total mixing
    }
    X.gaga = symmetric(X.gaga) # ensure contacts balance
    # print(X.gaga / aperm(X.gaga,c(3,4,1,2))) # DEBUG == 1 (exact) or NaN from 0/0
    # print(sweep(P.ga,2,rowSums(C.aa.y[[y]]),'*') / a.sum(X.gaga,c(3,4))) # DEBUG == 1 (approx)
    X.gaga.y[[y]] = X.gaga # copy result into output
  }
  names(X.gaga.y) = names(config$c.type)
  return(X.gaga.y) # "CX"
}

CX.norm = function(CX.gaga.y,P.ga,C.a.y=NULL){
  # convert CX (total contacts) to:
  # - Ci: per person, divide by P.ga (if C.ga.y is NULL)
  # - Cp: probability, divide by P.ga & C.ga (if C.ga.y is given)
  C.gaga.y = lapply(CX.gaga.y,function(X){ X*0 }) # initialize output
  for (y in seq(config$N$y)){
    for (g in seq(config$N$g)){
      for (a in seq(config$N$a)){
        c.a.y = ifelse(is.null(C.a.y),1,C.a.y[[y]][a])
        C.gaga.y[[y]][g,a,,] = CX.gaga.y[[y]][g,a,,] / P.ga[g,a] / c.a.y
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

main.mixing = function(t='REF',do.plot=TRUE){
  pop    = load.fsa.pop()
  C.AA.y = load.contacts()
  names(C.AA.y) = names(config$c.type)
  B.gg.t = load.group.mob(pop)
  C.aa.y = CAAy.to.Caay(C.AA.y)
  P.ga   = pop.to.Pga(pop)
  # TODO: apply g scaling
  CX.gaga.y = gen.mix.main(C.aa.y,P.ga,B.gg.t,t)
  Ci.gaga.y = CX.norm(CX.gaga.y,P.ga)
  Cp.gaga.y = CX.norm(CX.gaga.y,P.ga,lapply(C.aa.y,rowSums)) # unused for now
  if (do.plot){
    plot.contact.margins(C.AA.y,C.aa.y); ggsave(figname('C-restrat','contacts'),w=8,h=4)
    C.aa.y.diff = mapply(function(C.,C){ aggr.mix(C.,'Ci','a',P.ga) - C }, Ci.gaga.y, C.aa.y,SIMPLIFY=FALSE)
    plot.mix(C.aa.y,'Ci','a',trans='sqrt',clim=c(0,8.5));                        ggsave(figname('Caay','contacts'),w=8,h=4)
    plot.mix(C.aa.y.diff,'Ci','a',trans='nsqrt',clim=c(-2,+2),cmap='RdBu');      ggsave(figname('Caay-diff','contacts'),w=8,h=4)
    plot.mix(CX.gaga.y,'CX','a',P.ga=P.ga,trans='sqrt',aggr=TRUE);               ggsave(figname('CXaay','mixing',config$mode,t),w=8,h=4)
    plot.mix(CX.gaga.y,'CX','g',P.ga=P.ga,trans='sqrt',aggr=TRUE);               ggsave(figname('CXggy','mixing',config$mode,t),w=8,h=4)
    plot.mix(Ci.gaga.y,'Ci','a',P.ga=P.ga,trans='sqrt',aggr=TRUE,clim=c(0,8.5)); ggsave(figname('Ciaay','mixing',config$mode,t),w=8,h=4)
    plot.mix(Ci.gaga.y,'Ci','g',P.ga=P.ga,trans='sqrt',aggr=TRUE,clim=c(0,8.5)); ggsave(figname('Ciggy','mixing',config$mode,t),w=8,h=4)
    # plot.mix(Cp.gaga.y,'Cp','a',P.ga=P.ga,trans='sqrt'); ggsave(figname('Cpaay','mixing',config$mode,t),width=8,height=4)
    # plot.mix(Cp.gaga.y,'Cp','g',P.ga=P.ga,trans='sqrt'); ggsave(figname('Cpggy','mixing',config$mode,t),width=8,height=4)
  }
  return(Ci.gaga.y)
}
