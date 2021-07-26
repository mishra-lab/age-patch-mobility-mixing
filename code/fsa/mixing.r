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
  B.g.0 = rowSums(B.gg.t[['REF']]) # total mobility during REF
  B.g.t = rowSums(B.gg.t[[t]])     # total mobility during t
  B.gg = B.gg.t[[t]] + (B.g.t / B.g.0) * diag(1-B.g.0) # add mobility in same "g"
  B.hh = diag(N$g) # dummy "mobillity" for household contacts
  a.size = bin.size(config$age,norm=TRUE) # expected age distribution
  X.gaga.y = list() # initialize output
  for (y in seq(N$y)){ # for each contact type
    X.gaga = array(0,c(N$g,N$a,N$g,N$a),dimnames=config$X.names) # initialize X for this y
    C.gaga = aperm(replicate(N$g,replicate(N$g,C.aa.y[[y]])),c(3,1,4,2)) # pad C.aa to C.gaga
    # B.gg. = (1-h.y[y]) * B.gg + h.y[y] * B.hh # DEBUG: old method - similar to Arenas 2020
    for (g in seq(N$g)){ # for each group
      # traveller pool
      # Pt = P.ga * B.gg.[,g] * RC.g.y[[y]] # DEBUG: old method - similar to Arenas 2020
      Pt  = P.ga * (1-h.y[y]) * B.gg[,g] * RC.g.y[[y]] # pop of each "g" mixing in this pool
      PRt = sweep(Pt,2,a.size,'/')/sum(Pt) # normalize overall & age as C.gaga already weighted by age
      PRt[is.na(PRt)] = 0
      X = outer(Pt,PRt) * C.gaga # proportionate mixing * age mixing
      X.gaga = X.gaga + X # add contribution of this pool to the total mixing
      # home pool
      Ph  = P.ga[g,] * h.y[y] * RC.g.y[[y]] # pop (age groups) of this "g" at home
      PRh = Ph/a.size/sum(Ph) # normalize overall & age as C.aa.y already weighted by age
      PRh[is.na(PRh)] = 0
      X = outer(Ph,PRh) * C.aa.y[[y]] # proportionate mixing * age mixing
      X.gaga[g,,g,] = X.gaga[g,,g,] + X # add contribution of this pool to the total mixing
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

main.mixing = function(figdir='mx'){
  clim=c(0,12.3)
  # load stuff
  config   = set.config(mode='10x10',n.y='4')
  load(root.path('data','fsa','.rdata','Prem2017.rdata')) # -> Prem2017
  pop      = load.fsa.pop()
  P.ga     = pop.to.Pga(pop)
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
  P.A = data.frame(y=Prem2017$P.a,x=midpoint(config$age.contact))
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
  plot.mix(C.AA.y.0,'Ci','a',trans='sqrt',clim=clim); ggsave(figname('C4AAy0',figdir),w=14,h=4)
  plot.mix(C.AA.y.1,'Ci','a',trans='sqrt',clim=clim); ggsave(figname('C4AAy1',figdir),w=14,h=4)
  plot.mix(C.AA.y,  'Ci','a',trans='sqrt',clim=clim); ggsave(figname('C4AAy', figdir),w=14,h=4)
  plot.mix(C.aa.y,  'Ci','a',trans='sqrt'); ggsave(figname('C4aay', figdir),w=14,h=4)
  plot.mix(C.11.y,  'Ci','a',trans='sqrt') +
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
  B.gg.t = load.group.mob(pop)
  B.g.0 = rowSums(B.gg.t[['REF']])
  B.ggd.t = lapply(B.gg.t,function(B.gg){
    B.gg + pmin(1, rowSums(B.gg) / B.g.0) * diag(1-B.g.0)
  })
  plot.mix(B.gg.t[['REF']],'B','g',trans='sqrt');  ggsave(figname('Bgg',  'mx'),w= 5,h=4)
  plot.mix(B.ggd.t[['REF']],'B','g',trans='sqrt'); ggsave(figname('Bggd', 'mx'),w= 5,h=4)
  plot.mix(B.gg.t,'B','g',trans='sqrt');           ggsave(figname('Bggt', 'mx'),w=11,h=4)
  plot.mix(B.ggd.t,'B','g',trans='sqrt');          ggsave(figname('Bggdt','mx'),w=11,h=4)
  # compute 4D
  RC.g.y = gen.RC.g.y()
  CX.gaga.y = gen.mix.main(P.ga,C.aa.y,RC.g.y,B.gg.t,'REF')
  Ci.gaga.y = CX.norm(CX.gaga.y,P.ga)
  Ci.gaga   = Reduce('+',Ci.gaga.y)
  # plot 4D margins in a/a', g/g'
  plot.mix(Ci.gaga.y,'Ci','a',P.ga=P.ga,trans='sqrt',aggr=TRUE)
    ggsave(figname('CX4aay',figdir),w=14,h=4)
  plot.mix(Ci.gaga.y,'Ci','g',P.ga=P.ga,trans='sqrt',aggr=TRUE)
    ggsave(figname('CX4ggy',figdir),w=14,h=4)
  plot.mix(Ci.gaga,'Ci','a',P.ga=P.ga,trans='sqrt',aggr=TRUE)
    ggsave(figname('CXaay',figdir),w=5,h=4)
  plot.mix(Ci.gaga,'Ci','g',P.ga=P.ga,trans='sqrt',aggr=TRUE)
    ggsave(figname('CXggy',figdir),w=5,h=4)
  # compare 2D to 4D: should be exactly the same for REF
  P.a = replicate(config$N$a,colSums(P.ga)/mean(P.ga))
  C.aa.y. = lapply(C.aa.y,'*',P.a)
  C.aa.y.diff = mapply(function(C4,C2){ aggr.mix(C4,'Ci','a',P.ga) - C2 }, Ci.gaga.y, C.aa.y., SIMPLIFY=FALSE)
  plot.mix(C.aa.y.diff,'Ci','a',trans='nsqrt',clim=c(-1,+1),cmap='RdBu');
    ggsave(figname('C4aay-diff',figdir),w=14,h=4) # should be white
  # void plots for tikz figure (grabstract)
  load(root.path('data','fsa','.rdata','Prem2017.rdata'))
  C.list = list('0'=C.AA.y.0,'1'=C.AA.y.1,'2'=C.AA.y)
  P.list = list('0'=Prem2017$P.a/mean(Prem2017$P.a),'1'=rep(1,16),'2'=rep(1,16))
  lapply(names(C.list),function(n1){
    lapply(names(C.list[[n1]]),function(n2){
      C = C.list[[n1]][[n2]]
      plot.mix(C,'Ci','a',trans='sqrt') + void() + guides(fill='none')
        ggsave(root.path('code','tikz','mx','grabs','fig',paste0('CAA-',n1,'-',n2,'.pdf')),w=1,h=1)
    })
  })
  plot.mix(B.gg.t[['REF']],'B','g',trans='sqrt') + void() + guides(fill='none')
    ggsave(root.path('code','tikz','mx','grabs','fig','Bgg.pdf'),w=1,h=1)
  return(NULL)
}
