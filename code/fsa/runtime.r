fsa.mix.data = '.fsa-mix.rdata'

gen.mix.runtime.data = function(){
  source('data.r')
  source('mixing.r')
  pop    = load.fsa.pop()
  C.aa.y = load.contacts()
  B.gg   = load.group.mob(pop)
  C.a.y  = Caay.resample(C.aa.y)
  P.ga   = pop.to.Pga(pop)
  B.gg   = Reduce('+',B.gg[mo.ref]) / length(mo.ref)
  save(B.gg,P.ga,C.a.y,eps.y,h.y,file=fsa.mix.data)
}

a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}

gen.Ci.gg.y = function(f.B=1, c.mean=c(1,1), c.scale=c(.15,.30)){
  load(fsa.mix.data)
  N = list(y=2,g=10,a=12) # HACK
  B.tgg = f.B*B.gg + diag(1-rowSums(B.gg))
  Q.ga.y = lapply(seq(N$y),function(y){
    return(P.ga * outer(c.mean[y] + seq(+c.scale[y],-c.scale[y],l=10), C.a.y[[y]]))
  })
  X.gaga.y = list()
  for (y in seq(N$y)){
    X.y = array(0,c(N$g,N$a,N$g,N$a))
    for (x in seq(N$g)){ # travel-related
      Q.ga = Q.ga.y[[y]] * B.tgg[,x]
      X.yx = outer(Q.ga,Q.ga)/sum(Q.ga)
      Xa   = a.sum(X.yx,4)
      X.yx = X.yx * (1-eps.y[[y]])
      for (a in seq(N$a)){
        X.yx[,a,,a] = X.yx[,a,,a] + eps.y[[y]] * Xa[,a,]
      }
      X.y = X.y + (1-h.y[[y]]) * X.yx
    }
    for (g in seq(N$g)){ # non-travel-related
      Q.ga = Q.ga.y[[y]][g,]
      X.yg = outer(Q.ga,Q.ga) / sum(Q.ga)
      Xa   = a.sum(X.yg,2)
      X.yg = X.yg * (1-eps.y[[y]])
      for (a in seq(N$a)){
        X.yg[a,a] = X.yg[a,a] + eps.y[[y]] * Xa[a]
      }
      X.y[g,,g,] = X.y[g,,g,] + (h.y[[y]]) * X.yg
    }
    for (g in seq(N$g)){
      for (a in seq(N$a)){
        X.y[g,a,,] = X.y[g,a,,] / P.ga[g,a]
      }
    }
    X.gaga.y[[y]] = X.y
  }
  return(X.gaga.y)
}