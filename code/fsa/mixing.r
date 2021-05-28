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

Caay.to.Cay = function(C.y,age=NULL){
  if (is.null(age)){ age = info$age }
  mp = function(a){ a+c(diff(a)/2,0) }
  return(do.call(cbind,lapply(C.y,function(c.y){
    c.y.m = rowSums(c.y)
    return(approx(x=mp(age.contact),y=c.y.m,xo=mp(age),
      method='linear',yleft=c.y.m[1],yright=c.y.m[length(c.y.m)])$y)
  })))
}

gaga.mix = function(C.ay,P.ga,B.gg,mo){
  B.gg[['ref']] = Reduce('+',B.gg[mo.ref]) / length(mo.ref)
  B.tgg = B.gg[[mo]] + diag(1-rowSums(B.gg[['ref']])) # n.b. does not sum to one unless mo = ref
  Q.ga.y = lapply(info$c.type,function(y){ sweep(P.ga,2,C.ay[,y],'*') })
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

pct.self = function(X){
  return(median(diag(X) / rowSums(X)))
}

plot.mixing = function(B.gg,X.gaga,X.gaga.y,t='ref',f=NULL){
  if (is.null(f)){ f = file.path('mixing',MODE) }
  X.gaga = Reduce('+',X.gaga.y)
  print(pct.self(a.sum(X.gaga,c(2,4))))
  print(sapply(X.gaga.y,function(X){ pct.self(a.sum(X,c(2,4))) }))
  B.gg = append(B.gg,list('REF'=Reduce('+',B.gg[mo.ref]) / length(mo.ref)),1)
  for (mo in mo.ref){ B.gg[[mo]] = NULL }
  plot.mix(B.gg,    'g',clim=c(0,NA),aggr=FALSE); ggsave(figname('mobility',f),  width=10,height=8)
  plot.mix(X.gaga,  'g',clim=c(0,NA));            ggsave(figname('Xgg',     f,t),width= 5,height=4)
  plot.mix(X.gaga,  'g',clim=c(0,NA),xfun=offd);  ggsave(figname('Xgg-o',   f,t),width= 5,height=4)
  plot.mix(X.gaga,  'a',clim=c(0,NA));            ggsave(figname('Xaa',     f,t),width= 5,height=4)
  plot.mix(X.gaga,  'a',clim=c(0,NA),xfun=offd);  ggsave(figname('Xaa-o',   f,t),width= 5,height=4)
  plot.mix(X.gaga.y,'g',clim=c(0,NA));            ggsave(figname('Xggy',    f,t),width= 8,height=4)
  plot.mix(X.gaga.y,'g',clim=c(0,NA),xfun=offd);  ggsave(figname('Xggy-o',  f,t),width= 8,height=4)
  plot.mix(X.gaga.y,'a',clim=c(0,NA));            ggsave(figname('Xaay',    f,t),width= 8,height=4)
  plot.mix(X.gaga.y,'a',clim=c(0,NA),xfun=offd);  ggsave(figname('Xaay-o',  f,t),width= 8,height=4)
}

mixing.fname = function(t,sub='.raw/mix'){
  return(root.path('data','fsa',sub,paste0('mix_',MODE,'_',t,'.csv')))
}

save.mixing = function(X.gaga.y,t){
  mix. = do.call(expand.grid,dimnames(X.gaga.y[[1]]))
  mix = do.call(rbind,lapply(names(X.gaga.y),function(y){
    cbind(mix.,y=y,X=as.vector(X.gaga.y[[y]]))
  }))
  colnames(mix) = c('index.decile','index.age','other.decile','other.age','contact.type','n.contact')
  write.csv(mix,mixing.fname(t),row.names=FALSE)
}

merge.save.mixing = function(){
  mix = do.call(rbind,lapply(c('ref',mo.covid),function(t){
    mix. = cbind(month=t,read.csv(mixing.fname(t)))
  }))
  write.csv(mix,mixing.fname('all',''),row.names=FALSE)
}

main.mixing = function(t='ref'){
  C.y  = load.contacts()
  pop  = load.fsa.pop()
  B.gg = load.group.mob(pop)
  C.ay = Caay.to.Cay(C.y)
  P.ga = pop.to.Pga(pop)
  X.gaga.y = gaga.mix(C.ay,P.ga,B.gg,t)
  plot.mixing(B.gg,X.gaga,X.gaga.y,t)
  save.mixing(X.gaga.y,t)
}