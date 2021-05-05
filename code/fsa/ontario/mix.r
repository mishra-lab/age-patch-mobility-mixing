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

gaga.mix = function(C.ay,P.ga,B.gg,h.y){
  B.gg  = 0.25 * B.gg / rowSums(B.gg)
  B.g   = rowSums(B.gg)
  B.gg1 = B.gg + diag(1-B.g)
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