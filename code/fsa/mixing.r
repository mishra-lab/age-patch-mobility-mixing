source('data.r')
source('plot.r')

a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}

pop.to.Pga = function(pop){
  P = aggregate(pop~group+age,pop,sum)
  P.ga = matrix(P[order(P$age,P$group),]$pop,nrow=N$g,
    dimnames=list(g=names(info$group),a=names(info$age)))
  return(P.ga)
}

pmc.to.Cay = function(pmc){
  return(do.call(cbind,lapply(pmc,rowSums)))
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

plot.mixing = function(B.gg,X.gaga,X.gaga.y,pmc,t='ref',f='mixing'){
  X.gaga = Reduce('+',X.gaga.y)
  print(pct.self(a.sum(X.gaga,c(2,4))))
  print(sapply(X.gaga.y,function(X){ pct.self(a.sum(X,c(2,4))) }))
  B.gg = append(B.gg,list('REF'=Reduce('+',B.gg[mo.ref]) / length(mo.ref)),1)
  for (mo in mo.ref){ B.gg[[mo]] = NULL }
  # TODO: need clim per mode
  plot.mix(B.gg,    'g',clim=c(0,.21),aggr=FALSE); ggsave(figname('mobility',f),  width=10,height=8)
  plot.mix(X.gaga,  'g',clim=c(0, 12));            ggsave(figname('Xgg',     f,t),width= 5,height=4)
  plot.mix(X.gaga,  'g',clim=c(0,2.5),xfun=offd);  ggsave(figname('Xgg-o',   f,t),width= 5,height=4)
  plot.mix(X.gaga,  'a',clim=c(0, 35));            ggsave(figname('Xaa',     f,t),width= 5,height=4)
  plot.mix(X.gaga,  'a',clim=c(0, 15),xfun=offd);  ggsave(figname('Xaa-o',   f,t),width= 5,height=4)
  plot.mix(X.gaga.y,'g',clim=c(0,  8));            ggsave(figname('Xggy',    f,t),width= 8,height=4)
  plot.mix(X.gaga.y,'g',clim=c(0,2.5),xfun=offd);  ggsave(figname('Xggy-o',  f,t),width= 8,height=4)
  plot.mix(X.gaga.y,'a',clim=c(0, 30));            ggsave(figname('Xaay',    f,t),width= 8,height=4)
  plot.mix(pmc,     'a',clim=c(0, 10),aggr=FALSE); ggsave(figname('polymod', f,t),width= 8,height=4)
}

mixing.fname = function(t,sub='mix'){
  return(root.path('data','fsa',sub,paste0('mix_',mode,'x',mode,'_',t,'.csv')))
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
  pmc  = load.polymod()
  pop  = load.fsa.pop()
  B.gg = load.group.mob(pop)
  C.ay = pmc.to.Cay(pmc)
  P.ga = pop.to.Pga(pop)
  X.gaga.y = gaga.mix(C.ay,P.ga,B.gg,t)
  # plot.mixing(B.gg,X.gaga,X.gaga.y,pmc,t)
  save.mixing(X.gaga.y,t)
}