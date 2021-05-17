source('data.r')
source('plot.r')

a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}

pop.to.Pga = function(pop){
  G.n  = matrix(pop$decile,nrow=length(info$age))[1,] # WARNING: not robust
  S.gn = outer(info$decile,G.n,'==')*1
  P.na = matrix(pop$pop_agefsa,nrow=length(info$age)) # WARNING: not robust
  P.ga = S.gn %*% t(P.na)
  colnames(P.ga) = names(info$age)
  return(P.ga)
}

pmc.to.Cay = function(pmc){
  return(do.call(cbind,lapply(pmc,rowSums)))
}

gaga.mix = function(C.ay,P.ga,B.gg,mo){
  B.gg[['ref']] = Reduce('+',B.gg[mo.ref]) / length(mo.ref)
  B.gg1 = B.gg[[mo]] + diag(1-rowSums(B.gg[['ref']])) # n.b. does not sum to one unless B.gg = B.gg.ref
  Q.ga.y = lapply(info$c.type,function(y){ sweep(P.ga,2,C.ay[,y],'*') })
  X.gaga.y = list()
  for (y in seq(N$y)){
    X.y = array(0,c(N$g,N$a,N$g,N$a),dimnames=X.names)
    for (t in info$decile){
      Q.ga = Q.ga.y[[y]] * B.gg1[,t]
      T = sum(Q.ga)
      X.yt = outer(Q.ga,Q.ga)/T
      Xa   = a.sum(X.yt,4)
      X.yt = X.yt * (1-eps[[y]])
      for (a in seq(N$a)){ X.yt[,a,,a] = X.yt[,a,,a] + eps[[y]] * Xa[,a,] }
      X.y = X.y + X.yt
    }
    Xg = a.sum(X.y,3)
    X.y = X.y * (1-h.y[[y]])
    for (g in seq(N$g)){ X.y[g,,g,] = X.y[g,,g,] + h.y[[y]] * Xg[g,,] }
    X.gaga.y[[y]] = X.y
  }
  names(X.gaga.y) = names(info$c.type)
  # y = 6; print(a.sum(X.gaga.y[[y]],c(3,4))/Q.ga.y[[y]]) # DEBUG
  return(X.gaga.y)
}

pct.self = function(X){
  return(median(diag(X) / rowSums(X)))
}

main.mixing = function(t='ref'){
  pmc  = load.polymod()
  pop  = load.fsa.pop()
  B.gg  = load.decile.mob()
  C.ay = pmc.to.Cay(pmc)
  P.ga = pop.to.Pga(pop)
  X.gaga.y = gaga.mix(C.ay,P.ga,B.gg,t)

  X.gaga = Reduce('+',X.gaga.y)
  print(pct.self(a.sum(X.gaga,c(2,4))))
  print(sapply(X.gaga.y,function(X){ pct.self(a.sum(X,c(2,4))) }))

  g = plot.mix(B.gg,    'g',clim=c(0,.15),aggr=FALSE); ggsave(figname('mobility'),width= 8,height= 8)
  g = plot.mix(X.gaga,  'g',clim=c(0, 18));            ggsave(figname('Xgg'),     width= 5,height= 4)
  g = plot.mix(X.gaga,  'g',clim=c(0,  1),xfun=offd);  ggsave(figname('Xgg-o'),   width= 5,height= 4)
  g = plot.mix(X.gaga,  'a',clim=c(0, 32));            ggsave(figname('Xaa'),     width= 5,height= 4)
  g = plot.mix(X.gaga,  'a',clim=c(0, 13),xfun=offd);  ggsave(figname('Xaa-o'),   width= 5,height= 4)
  g = plot.mix(X.gaga.y,'g',clim=c(0,  6));            ggsave(figname('Xggy'),    width=14,height= 3)
  g = plot.mix(X.gaga.y,'g',clim=c(0,.50),xfun=offd);  ggsave(figname('Xggy-o'),  width=14,height= 3)
  g = plot.mix(X.gaga.y,'a',clim=c(0, 10));            ggsave(figname('Xaay'),    width=14,height= 3)
  g = plot.mix(pmc,     'a',clim=c(0,  5),aggr=FALSE); ggsave(figname('polymod'), width=14,height= 3)
}