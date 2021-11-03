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

decile.pop.to.Pga = function(){
  # pop: dataframe with columns (at least): group, age, pop(ulation)
  # P.ga: [value] matrix of population by group & age
  P = read.csv(root.path('data', 'pop_decile.csv'))
  P.ga = matrix(P[order(P$age,P$group),]$pop,nrow=config$N$g,
                dimnames=list(g=names(config$group),a=names(config$age)))
  return(P.ga)
}

resample.contacts = function(Ci,xi,xo){
  # Ci: contact matrix stratified by xi & xi (square)
  # xi: list of cut points in Ci (lower cut points of each strata)
  # xo: list of cut points in Co (lower cut points of each strata)
  # Co: [value] contact matrix re-stratified by xo * xo (square)
  # We cannot simply interp2d, since the number of contacts per person would change;
  # that's why we interpolate to m1 * m1 first, then aggregate using A
  # and normalize by the strata sizes (bin.size);
  # the number of contacts per person might change a little;
  # we also pad Ci diagonally before interp to avoid any NA or edge effects.
  xie = bin.extrap(xi,x.max=max(bin.extrap(xo)))
  xip = bin.extrap(xie,pre=1,post=1)
  Cip = diag.pad(Ci,pre=1,post=1+length(xie)-length(xi))
  x1  = seq(min(xip),max(xip))
  m1  = midpoint(x1)
  C1  = interp2d(Cip/bin.size(xip),midpoint(xip),m1)
  i   = (m1 > xie[1]) & (x1 < xie[length(xie)])
  # plot.mix(C1[i,i],'Ci','a'); ggsave('Rplots.pdf'); # DEBUG
  A  = tail(outer(c(0,xo),m1[i],'<') * outer(c(xo,Inf),m1[i],'>'),-1)
  Co = (A %*% C1[i,i] %*% t(A)) / bin.size(xo)
  return(Co)
}

CAAy.to.Caay = function(C.AA.y,age.out=NULL){
  # C.AA.y: list of contact matrices in original age stratification
  # age.out: list of cut points in C.aa.y (desired age stratification)
  # C.aa.y: [value] list of contact matrices in desired age stratification
  if (is.null(age.out)){ age.out = config$age }
  return(lapply(C.AA.y,resample.contacts,xi=config$age.contact,xo=age.out))
}

gen.RC.g.y = function(){
  # RC.g.y: [value] list by contact type of vectors: relative C by group
  return(lapply(config$c.type,function(c.type){
    return(config$RC.global[c.type] * config$RC.decile/mean(config$RC.decile))
  }))
}

gen.mix.main = function(P.ga,C.aa.y,RC.g.y,B.gg.t,t){
  # estimates the total number of contacts between all combinations of age/geo groups
  # C.ga.y: list by contact type of matrices: contacts by group & age
  # P.ga:   matrix of population by group & age
  # B.gg.t: list by time of matrices: mobility by group & other group
  # t:      month
  # X.gaga.y: [value] list by type of arrays: total number of contacts
  #           by age, group, other age, other group (not per person)
  N = config$N; h.y = config$h.y; # convenience
  B.gg = B.gg.t[[t]] # mobility during t
  B.hh = diag(N$g) # dummy "mobillity" for household contacts
  a.size = bin.size(config$age,norm=TRUE) # expected age distribution
  if (config$method$age.adapt){
    P.ga. = P.ga
  } else {
    P.ga. = replicate(N$a,rowMeans(P.ga))
    dimnames(P.ga.) = dimnames(P.ga)
  }
  if (! config$method$age.by.type){
    C.aa.u = Reduce('+',C.aa.y)
    C.aa.y = lapply(C.aa.y,function(C){
      C.aa.u * rowSums(C)/rowSums(C.aa.u)
    })
  }
  # PR = replicate(N$a,rowMeans(sweep(P.ga,2,a.size,'/')/sum(P.ga)))
  X.gaga.y = list() # initialize output
  for (y in seq(N$y)){ # for each contact type
    X.gaga = array(0,c(N$g,N$a,N$g,N$a),dimnames=config$X.names) # initialize X for this y
    C.gaga = aperm(replicate(N$g,replicate(N$g,C.aa.y[[y]])),c(3,1,4,2)) # pad C.aa to C.gaga
    for (g in seq(N$g)){ # for each group
      # traveller pool
      if (config$method$home.pool){ # pop of each "g" mixing in this pool ...
        Pt  = RC.g.y[[y]] * P.ga  * (1-h.y[y]) * B.gg[,g]
        Pt. = RC.g.y[[y]] * P.ga. * (1-h.y[y]) * B.gg[,g]
      } else {
        Pt  = RC.g.y[[y]] * P.ga  * ((1-h.y[y]) * B.gg[,g] + (h.y[y]) * B.hh[,g])
        Pt. = RC.g.y[[y]] * P.ga. * ((1-h.y[y]) * B.gg[,g] + (h.y[y]) * B.hh[,g])
      }
      pt. = sweep(Pt.,2,a.size,'/')/sum(Pt.) # normalize overall & age as C.gaga already weighted by age
      pt.[is.nan(pt.)] = 0
      X = outer(Pt,pt.) * C.gaga # proportionate mixing * age mixing
      X.gaga = X.gaga + X # add contribution of this pool to the total mixing
      # home pool
      if (config$method$home.pool){
        Ph  = RC.g.y[[y]] * P.ga[g,]  * h.y[y] # pop (age groups) of this "g" at home
        Ph. = RC.g.y[[y]] * P.ga.[g,] * h.y[y]
        ph. = Ph./a.size/sum(Ph)# normalize overall & age as C.aa.y already weighted by age
        ph.[is.nan(ph.)] = 0
        X = outer(Ph,ph.) * C.aa.y[[y]] # proportionate mixing * age mixing
        X.gaga[g,,g,] = X.gaga[g,,g,] + X # add contribution of this pool to the total mixing
      }
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
  for (y in seq(length(CX.gaga.y))){
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

melt.mixing = function(X.gaga.y.t,what){
  # TODO: use what
  X. = do.call(expand.grid,dimnames(X.gaga.y.t[[1]][[1]]))
  return(do.call(rbind,lapply(names(X.gaga.y.t),function(t){
    do.call(rbind,lapply(seq(config$N$y),function(y){
      cbind(t=t,y=config$c.type[y],X.,C=as.vector(X.gaga.y.t[[t]][[y]]))
    }))
  })))
}

gen.save.mixing = function(do.save=TRUE,norm=TRUE,t=NULL,...){
  config = set.config(...)
  P.ga   = decile.pop.to.Pga()
  C.AA.y = load.contacts()
  C.aa.y = CAAy.to.Caay(C.AA.y)
  RC.g.y = gen.RC.g.y()
  B.gg.t = load.group.mob(pop,B='B')
  if (missing(t)){ t = names(B.gg.t) }
  Ci.gaga.y.t = lapply(t,function(ti){
    CX.gaga.y = gen.mix.main(P.ga,C.aa.y,RC.g.y,B.gg.t,ti)
    if (norm){
      CX.gaga.y = CX.norm(CX.gaga.y,P.ga)
    }
    return(CX.gaga.y)
  })
  names(Ci.gaga.y.t) = t
  if (do.save){
    X = melt.mixing(Ci.gaga.y.t)
    write.csv(X,root.path('data','Ci_gagayt.csv'),row.names=FALSE)
  } else {
    return(Ci.gaga.y.t)
  }
}

compare.mixing = function(figdir='compare'){
  set.config(n.y = '4',
             age = config$age.contact,
             h.y = c('home'=1,'work'=0.1,'school'=0.5,'others'=0.3))
  P.ga   = decile.pop.to.Pga()
  C.AA.y = load.contacts()
  C.aa.y = CAAy.to.Caay(C.AA.y)
  RC.g.y = gen.RC.g.y()
  B.gg.t = load.group.mob(pop)
  mix.fun = function(key,t='REF',...){
    args = list(...)
    for (name in names(args))
    config$method[[name]] <<- args[[name]]
    CX = CX.norm(gen.mix.main(P.ga, C.aa.y, RC.g.y, B.gg.t, t),P.ga)
    CX. = Reduce('+',CX)
    plot.mix(CX,'Ci','a',P.ga=P.ga,aggr=TRUE); ggsave(figname(paste0('C4iaa',key),figdir),w=14,h=4)
    plot.mix(CX,'Ci','g',P.ga=P.ga,aggr=TRUE); ggsave(figname(paste0('C4igg',key),figdir),w=14,h=4)
    plot.mix(CX.,'Ci','a',P.ga=P.ga,aggr=TRUE); ggsave(figname(paste0('Ciaa',key),figdir),w=5,h=4)
    plot.mix(CX.,'Ci','g',P.ga=P.ga,aggr=TRUE); ggsave(figname(paste0('Cigg',key),figdir),w=5,h=4)
    return(CX)
  }
  Ci.gaga.y.m = list( # steps: home pool, age adapt, age by type
    '3' = mix.fun('3',home.pool=TRUE, age.adapt=TRUE, age.by.type=TRUE),
    '2' = mix.fun('2',home.pool=FALSE,age.adapt=TRUE, age.by.type=TRUE),
    '1' = mix.fun('1',home.pool=FALSE,age.adapt=FALSE,age.by.type=TRUE),
    '0' = mix.fun('0',home.pool=FALSE,age.adapt=FALSE,age.by.type=FALSE))
  comb = combn(names(Ci.gaga.y.m),2)
  max.abs = function(x){ max(abs(x)) }
  cmfun = function(cm){ cm = max(cm,.2); return(c(-cm,+cm)) }
  for (m in c('a','g')){
    for (i in seq(ncol(comb))){
      Ci.y.diff = lapply(names(config$c.type),function(y){
        aggr.mix(Ci.gaga.y.m[[comb[1,i]]][[y]],'Ci',m,P.ga) -
        aggr.mix(Ci.gaga.y.m[[comb[2,i]]][[y]],'Ci',m,P.ga)
      })
      names(Ci.y.diff) = names(config$c.type)
      cm = max(sapply(Ci.y.diff,max.abs))
      plot.mix(Ci.y.diff,'Ci',m,P.ga=P.ga,cmap='RdBu',trans='nsqrt',gez=FALSE,clim=cmfun(cm))
        ggsave(figname(paste0('D4',m,m,comb[1,i],'v',comb[2,i]),figdir),w=14,h=4)
      Ci.diff = Reduce('+',Ci.y.diff)
      cm = max.abs(Ci.diff)
      plot.mix(Ci.diff,'Ci',m,P.ga=P.ga,cmap='RdBu',trans='nsqrt',gez=FALSE,clim=cmfun(cm))
        ggsave(figname(paste0('D',m,m,comb[1,i],'v',comb[2,i]),figdir),w=5,h=4)
    }
  }
}

main.mixing = function(figdir=''){
  clim=c(0,12.3)
  # load stuff
  config   = set.config(mode='10x10',n.y='4')
  load(root.path('data','.rdata','Prem2021.rdata'))
  P.ga   = decile.pop.to.Pga()
  P.a      = colSums(P.ga)
  C.AA.y   = load.contacts()
  C.AA.y.0 = load.contacts(P.norm=FALSE,sym=FALSE)
  C.AA.y.1 = load.contacts(P.norm=TRUE, sym=FALSE)
  # do restratify
  age.1  = seq(100); names(age.1) = age.1;
  C.11.y = CAAy.to.Caay(C.AA.y,age.1)
  C.aa.y = CAAy.to.Caay(C.AA.y)
  C.aa.y.p = lapply(C.aa.y,function(C){ sweep(C,2,P.a,'*') / weighted.mean(P.a,bin.size(config$age)) })
  # plot age distribution
  P.A = data.frame(y=Prem2021$P.a,x=midpoint(config$age.contact))
  P.a = data.frame(y=colMeans(P.ga),x=midpoint(config$age))
  plot.pop.density() +
    geom_line(data=P.A,aes(x=x,y=y/mean(y)),linetype='dashed')
  ggsave(figname('Pga',figdir),w=8,h=4)
  # compare margins
  C.list   = list(C.AA.y.0, C.AA.y.1, C.AA.y, C.aa.y, C.aa.y.p)
  age.list = list(config$age.contact,config$age.contact,config$age.contact,config$age,config$age)
  names(C.list) = c('Original','Unweighted','Unweighted\n+Balanced','Unweighted\nTarget','Target')
  names(age.list) = names(C.list)
  plot.contact.margins(C.list,age.list,style='dots') +
    guides(color=guide_legend(title='',keyheight=2)) +
    scale_x_continuous(minor_breaks=NULL,
      breaks=unname(bin.extrap(config$age.contact),minor),
      sec.axis=sec_axis(~.,breaks=unname(bin.extrap(config$age))))
    ggsave(figname('C4ay',figdir),w=14,h=4)
  # compare the various versions of C.AA.y
  b.1 = seq(0,length(age.1),5)
  plot.mix(C.AA.y.0,'Ci','a',clim=clim); ggsave(figname('C4AAy0',figdir),w=14,h=4)
  plot.mix(C.AA.y.1,'Ci','a',clim=clim); ggsave(figname('C4AAy1',figdir),w=14,h=4)
  plot.mix(C.AA.y,  'Ci','a',clim=clim); ggsave(figname('C4AAy', figdir),w=14,h=4)
  plot.mix(C.aa.y,  'Ci','a'); ggsave(figname('C4aay', figdir),w=14,h=4)
  plot.mix(C.11.y,  'Ci','a') +
    scale_x_discrete(breaks=b.1) +
    scale_y_discrete(breaks=b.1); ggsave(figname('C411y',figdir),w=14,h=4)
  # difference plots
  C.AA.y.d01 = mapply('-', C.AA.y.1, C.AA.y.0, SIMPLIFY=FALSE)
  C.AA.y.d12 = mapply('-', C.AA.y,   C.AA.y.1, SIMPLIFY=FALSE)
  C.AA.y.d02 = mapply('-', C.AA.y,   C.AA.y.0, SIMPLIFY=FALSE)
  plot.mix(C.AA.y.d01,'Ci','a',trans='nsqrt',clim=c(-2,+2),cmap='RdBu',gez=FALSE)
    ggsave(figname('C4AAy-d01',figdir),w=14,h=4)
  plot.mix(C.AA.y.d12,'Ci','a',trans='nsqrt',clim=c(-2,+2),cmap='RdBu',gez=FALSE)
    ggsave(figname('C4AAy-d12',figdir),w=14,h=4)
  plot.mix(C.AA.y.d02,'Ci','a',trans='nsqrt',clim=c(-2,+2),cmap='RdBu',gez=FALSE)
    ggsave(figname('C4AAy-d02',figdir),w=14,h=4)
  # mobility plots
  Bc.gg.t = load.group.mob(pop,B='Bc')
  Ba.gg.t = load.group.mob(pop,B='Ba')
  B.gg.t  = load.group.mob(pop,B='B')
  plot.mix(Bc.gg.t[['REF']],'Bc','g'); ggsave(figname('Bcgg', figdir),w= 5,h=4)
  plot.mix(Ba.gg.t[['REF']],'B','g');  ggsave(figname('Bagg', figdir),w= 5,h=4)
  plot.mix(B.gg.t[['REF']], 'B','g');  ggsave(figname('Bgg',  figdir),w= 5,h=4)
  plot.mix(Bc.gg.t,'Bc','g');          ggsave(figname('Bcggt',figdir),w=10,h=4)
  plot.mix(Ba.gg.t,'B','g');           ggsave(figname('Baggt',figdir),w=10,h=4)
  plot.mix(B.gg.t, 'B','g');           ggsave(figname('Bggt', figdir),w=10,h=4)
  plot.mix(Ba.gg.t,'B','g',trans='log10'); ggsave(figname('Baggt-log',figdir),w=10,h=4)
  plot.mix(B.gg.t, 'B','g',trans='log10'); ggsave(figname('Bggt-log', figdir),w=10,h=4)
  # compute 4D
  RC.g.y = gen.RC.g.y()
  CX.gaga.y = gen.mix.main(P.ga,C.aa.y,RC.g.y,B.gg.t,'REF')
  Ci.gaga.y = CX.norm(CX.gaga.y,P.ga)
  Ci.gaga   = Reduce('+',Ci.gaga.y)
  # plot 4D margins in a/a', g/g'
  plot.mix(Ci.gaga.y,'Ci','a',P.ga=P.ga,aggr=TRUE)
    ggsave(figname('CX4aay',figdir),w=14,h=4)
  plot.mix(Ci.gaga.y,'Ci','g',P.ga=P.ga,aggr=TRUE)
    ggsave(figname('CX4ggy',figdir),w=14,h=4)
  plot.mix(Ci.gaga,'Ci','a',P.ga=P.ga,aggr=TRUE)
    ggsave(figname('CXaay',figdir),w=5,h=4)
  plot.mix(Ci.gaga,'Ci','g',P.ga=P.ga,aggr=TRUE)
    ggsave(figname('CXggy',figdir),w=5,h=4)
  # compare 2D to 4D: should be exactly the same for REF
  P.a = replicate(config$N$a,colSums(P.ga)/mean(P.ga))
  C.aa.y. = lapply(C.aa.y,'*',P.a)
  C.aa.y.diff = mapply(function(C4,C2){ aggr.mix(C4,'Ci','a',P.ga) - C2 }, Ci.gaga.y, C.aa.y., SIMPLIFY=FALSE)
  plot.mix(C.aa.y.diff,'Ci','a',trans='nsqrt',clim=c(-1,+1),cmap='RdBu');
    ggsave(figname('C4aay-diff',figdir),w=14,h=4) # should be white
  # void plots for tikz figure (grabstract)
  tikzpath = function(...){ root.path('code','tikz','methodsx','grabstract','fig',...) }
  load(root.path('data','.rdata','Prem2021.rdata'))
  C.list = list('0'=C.AA.y.0,'1'=C.AA.y.1,'2'=C.AA.y)
  P.list = list('0'=Prem2021$P.a/mean(Prem2021$P.a),'1'=rep(1,16),'2'=rep(1,16))
  lapply(names(C.list),function(n1){
    lapply(names(C.list[[n1]]),function(n2){
      C = C.list[[n1]][[n2]]
      plot.mix(C,'Ci','a') + void() + guides(fill='none')
        ggsave(tikzpath(paste0('CAA-',n1,'-',n2,'.pdf')),w=1,h=1)
    })
  })
  plot.mix(Bc.gg.t[['REF']],'Bc','g') + void() + guides(fill='none')
    ggsave(tikzpath('Bcgg.pdf'),w=1,h=1)
  return(NULL)
}
