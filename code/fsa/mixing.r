source('utils.r')
source('data.r')
source('plot.r')

resample.contacts = function(C5,ai,ao){
  # WARNING: assumes C5 is 5-year age groups
  midpoints = function(a){ a+c(diff(a)/2,.5) }
  m1 = midpoints(seq(0,79))
  C1 = interp2d(C5,midpoints(ai),m1)/5
  M = t(tail(outer(c(0,ao),m1,'<') * outer(c(ao,Inf),m1,'>'),-1))
  Co = t(apply(C1,1,interp1d,xi=m1,xo=midpoints(ao)) %*% M)
  dimnames(Co) = list(a=names(ao),a.=names(ao))
  return(Co)
}

pop.to.Pga = function(pop){
  P = aggregate(pop~group+age,pop,sum)
  P.ga = matrix(P[order(P$age,P$group),]$pop,nrow=N$g,
    dimnames=list(g=names(info$group),a=names(info$age)))
  return(P.ga)
}

Caay.to.Cgay = function(C.aa.y){
  # WARNING: assume 10 deciles
  S = map.decile(data.frame(decile=seq(10)))
  mp = function(a){ a+c(diff(a)/2,0) }
  C.ga.y = lapply(names(C.aa.y),function(y){
    S$scale = f.c.mean[[y]]+seq(+f.c.scale[[y]],-f.c.scale[[y]],l=10)
    c.a.i = rowSums(C.aa.y[[y]])
    c.a.o = approx(x=mp(age.contact),y=c.a.i,xo=mp(info$age),
      method='linear',yleft=c.a.i[1],yright=c.a.i[length(c.a.i)])$y
    return(matrix(aggregate(CS~group+age,
        transform(merge(data.frame(C=c.a.o,age=info$age),S),CS=C*scale),mean)$CS,
      nrow=N$g,dimnames=X.names[1:2]))
  })
  names(C.ga.y) = names(C.aa.y)
  return(C.ga.y)
}

gen.mix.total = function(C.ga.y,P.ga,B.gg,mo){
  # estimates the total number of contacts between all combinations of age/geo groups
  # C.ga.y: contacts per group, age & contact type
  # P.ga: population per geo & age
  # B.gg: mobility between geo groups (home / other)
  # mo:   month (t)
  B.gg[['ref']] = Reduce('+',B.gg[mo.ref]) / length(mo.ref)
  B.tgg = B.gg[[mo]] + diag(1-rowSums(B.gg[['ref']])) # n.b. does not sum to one unless mo = ref
  Q.ga.y = lapply(info$c.type,function(y){ P.ga * C.ga.y[[y]] })
  X.gaga.y = list()
  for (y in seq(N$y)){
    X.y = array(0,c(N$g,N$a,N$g,N$a),dimnames=X.names)
    for (x in seq(N$g)){ # travel-related
      Q.ga = Q.ga.y[[y]] * B.tgg[,x]
      T    = sum(Q.ga)
      X.yx = outer(Q.ga,Q.ga)/T
      Xa   = a.sum(X.yx,4)
      X.yx = X.yx * (1-eps[[y]])
      for (a in seq(N$a)){ X.yx[,a,,a] = X.yx[,a,,a] + eps[[y]] * Xa[,a,] }
      X.y = X.y + (1-h.y[[y]]) * X.yx
    }
    for (g in seq(N$g)){ # non-travel-related
      Q.ga = Q.ga.y[[y]][g,]
      X.yg = outer(Q.ga,Q.ga) / sum(Q.ga)
      Xa   = a.sum(X.yg,2)
      X.yg = X.yg * (1-eps[[y]])
      for (a in seq(N$a)){ X.yg[a,a] = X.yg[a,a] + eps[[y]] * Xa[a] }
      X.y[g,,g,] = X.y[g,,g,] + (h.y[[y]]) * X.yg
    }
    X.gaga.y[[y]] = X.y
  }
  names(X.gaga.y) = names(info$c.type)
  # for (y in seq(N$y)){ print(a.sum(X.gaga.y[[y]],c(3,4))/Q.ga.y[[y]]) } # DEBUG
  return(X.gaga.y)
}

CX.norm = function(CX.gaga.y,P.ga,C.ga.y=NULL){
  C.gaga.y = lapply(CX.gaga.y,function(X){ X*0 })
  for (y in seq(N$y)){
    for (g in seq(N$g)){
      for (a in seq(N$a)){
        c.ga.y = ifelse(is.null(C.ga.y),1,C.ga.y[[y]][g,a])
        C.gaga.y[[y]][g,a,,] = CX.gaga.y[[y]][g,a,,] / P.ga[g,a] / c.ga.y
      }
    }
  }
  names(C.gaga.y) = names(CX.gaga.y)
  return(C.gaga.y)
}

aggr.mix = function(C.gaga,what,vs,P.ga=NULL,aggr=TRUE){
  if (!aggr){ return(C.gaga) }
  d = switch(vs,'a'=c(1,3),'g'=c(2,4),'i'=c(3,4)) # 'o'=c(1,2) TODO: check interpretation of 'o'
  if (what=='CX'){ return( a.sum(C.gaga,d) / 1e6 )  }
  if (what=='Ci'){
    P.gaga = array(P.ga,dim(C.gaga))
    return( a.mean( a.sum(C.gaga,d[d>=3]), d[d<=2], a.sum(P.gaga,d[d>=3])) )
  }
  if (what=='Cp'){ stop('aggr.mix for Cp not yet implemented') }
}

pct.self = function(X){
  return(median(diag(X) / rowSums(X)))
}

plot.mixing = function(C.gaga.y,what,t='ref',...){
  print(sapply(C.gaga.y,function(C){ pct.self( aggr.mix(C,what,'g',...) ) }))
  plot.mix(C.gaga.y,what,'g',...,clim=c(0,NA));            ggsave(figname(paste0(what,'ggy'),  'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'g',...,clim=c(0,NA),xfun=offd);  ggsave(figname(paste0(what,'ggy-o'),'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'a',...,clim=c(0,NA));            ggsave(figname(paste0(what,'aay'),  'mixing',MODE,t),width=8,height=4)
  plot.mix(C.gaga.y,what,'a',...,clim=c(0,NA),xfun=offd);  ggsave(figname(paste0(what,'aay-o'),'mixing',MODE,t),width=8,height=4)
}

mixing.fname = function(what,t,sub='.tmp/mix'){
  path = root.path('data','fsa',sub)
  suppressWarnings({ dir.create(path,recursive=TRUE) })
  return(file.path(path,paste0(what,'_',MODE,'_',t,'.csv')))
}

save.mixing = function(X.gaga.y,what,t){
  mix. = do.call(expand.grid,dimnames(X.gaga.y[[1]]))
  mix = do.call(rbind,lapply(names(X.gaga.y),function(y){
    cbind(mix.,y=info$c.type[[y]],X=as.vector(X.gaga.y[[y]]))
  }))
  v.name = switch(what,CX='n.contacts',Ci='n.contact.pp',Cp='prob.contact')
  colnames(mix) = c('index.decile','index.age','other.decile','other.age','contact.type',v.name)
  write.csv(mix,mixing.fname(what,t),row.names=FALSE)
}

merge.save.mixing = function(what){
  mix = do.call(rbind,lapply(c('ref',mo.covid),function(t){
    mix. = cbind(month=t,read.csv(mixing.fname(what,t)))
  }))
  write.csv(mix,mixing.fname(what,'all',''),row.names=FALSE)
}

main.mixing = function(t='ref'){
  pop    = load.fsa.pop()
  C.aa.y = load.contacts()
  B.gg   = load.group.mob(pop)
  C.ga.y = Caay.to.Cgay(C.aa.y)
  P.ga   = pop.to.Pga(pop)
  CX.gaga.y = gen.mix.total(C.ga.y,P.ga,B.gg,t) # absolute contacts
  Ci.gaga.y = CX.norm(CX.gaga.y,P.ga)           # contacts per person
  Cp.gaga.y = CX.norm(CX.gaga.y,P.ga,C.ga.y)    # contact probability
  plot.mixing(CX.gaga.y,'CX',t)
  plot.mixing(Ci.gaga.y,'Ci',t,P.ga=P.ga)
  save.mixing(Ci.gaga.y,'Ci',t)
}