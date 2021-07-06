suppressPackageStartupMessages({
library(reshape2)
library(viridis)
library(ggplot2)
library(ggridges)
library(scales)
})

# supporting functions for plotting stuff

offd  = function(x){ x-diag(rep(NA,nrow(x))) }

mix.melt = function(C,what,vs,aggr=TRUE,xfun=NULL,...){
  if (is.null(xfun)){ xfun = identity }
  if (is.list(C)){
    return(do.call(rbind,lapply(names(C),function(name){
      C.aggr = aggr.mix(C[[name]],what,vs,aggr=aggr,...)
      C. = melt(xfun(C.aggr),value.name='X',varnames=c('i','i.'))
      C.$group = factor(name,levels=names(C))
      return(C.)
    })))
  } else {
    C.aggr = aggr.mix(C,what,vs,aggr=aggr,...)
    return( melt(xfun(C.aggr),value.name='X',varnames=c('i','i.')) )
  }
}
nsqrt_trans = function(){
  trans_new('nsqrt',function(x){sign(x)*sqrt(abs(x))},function(x){sign(x)*x^2})
}
cmap.fun = function(do='color',cmap='inferno',discrete=FALSE,...){
  if (cmap %in% c('inferno')){
    return(get(paste('scale',do,'viridis',sep='_'))(
      option=cmap,end=.95,discrete=discrete,...))
  }
  if (cmap %in% c('Spectral','RdBu')){
    return(get(paste('scale',do,ifelse(discrete,'brewer','distiller'),sep='_'))(
      palette=cmap,...))
  }
}
plot.mix = function(C,what,vs,aggr=FALSE,xfun=NULL,trans='identity',clim=c(0,NA),cmap='inferno',...){
  C. = mix.melt(C,what,vs,aggr=aggr,xfun=xfun,...)
  v.name = switch(what,
    CX = 'Total Contacts\n(Millions)',
    Ci = 'Contacts\nPer Person',
    Cp = 'Contact\nFormation\nProbability',
    B  = '% Decile Pop\nwho Travelled\nto Other Decile')
  g = ggplot(C.,aes(x=factor(i.),y=factor(i),fill=X,color=X)) +
    geom_tile() +
    coord_fixed(ratio=1) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(x=config$labels[[vs]]$x,y=config$labels[[vs]]$y,fill=v.name) +
    cmap.fun('fill', cmap=cmap,limits=clim,na.value='transparent',trans=trans) +
    cmap.fun('color',cmap=cmap,limits=clim,na.value='transparent',trans=trans) +
    theme_light() +
    guides(color='none',fill=guide_colorbar(barheight=5)) +
    switch(vs,
      a = theme(axis.text.x=element_text(angle=90,hjust=1)),
      n = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    )
  if (is.list(C)){
    if (length(C) <= 6) { g = g + facet_grid(cols=vars(group)) }
    else { g = g + facet_wrap(vars(group),ncol=6) }
  }
  return(g)
}
plot.device.density = function(X,x,y='month',bw=NULL,q=4,xmax=NULL,legend=FALSE){
  if (is.null(xmax)){ xmax=max(X[[x]]) }
  if (y=='month'){ y.lab='Month' } else { y.lab = y }
  g = ggplot(X,aes_string(y=y,fill='factor(stat(quantile))',x=x)) +
    stat_density_ridges(
      geom='density_ridges_gradient',
      calc_ecdf=TRUE,
      quantiles=q,
      bandwidth=bw) +
    scale_fill_viridis(discrete=TRUE,option='inferno',alpha=.7,begin=.2,end=.8) +
    scale_y_discrete(limits=rev(levels(X[[y]]))) +
    xlim(0,xmax) + labs(x=x,y=y.lab,fill='Quantile') +
    theme_light()
  if (!legend){
    g = g + guides(fill='none')
  }
  return(g)
}
plot.pop.density = function(){
  pop = load.fsa.pop(aggr=FALSE)
  P = aggregate(pop~group+age,pop,sum)
  g = ggplot(P,aes(x=age,y=pop,group=group,color=group,fill=group)) +
    geom_density(stat='smooth',alpha=.1,span=.1) +
    cmap.fun('color','Spectral',discrete=TRUE) +
    cmap.fun('fill','Spectral',discrete=TRUE) +
    labs(x='Age',y='Population',color='Decile',fill='Decile') +
    theme_light()
  return(g)
}
plot.contact.margins = function(C.AA.y,C.aa.y){
  c.melt = function(C,age,vs){
    return(do.call('rbind',lapply(seq(2),function(y){
      return(data.frame(Contacts=rowSums(C[[y]]),Age=age,Type=names(config$c.type)[y],vs=vs))
    })))
  }
  C.A = c.melt(C.AA.y,age=midpoint(config$age.contact),vs='Original')
  C.a = c.melt(C.aa.y,age=midpoint(config$age),vs='Restratified')
  g = ggplot(map=aes(x=Age,y=Contacts,group=vs,color=vs)) +
    geom_line(data=C.A) + geom_point(data=C.A) +
    geom_line(data=C.a) + geom_point(data=C.a) +
    facet_grid(cols=vars(Type)) +
    theme_light()
  return(g)
}
plot.cases = function(cases){
  g = ggplot(cases,aes(x=as.Date(t),y=cases,group=decile,color=decile)) +
    geom_line() +
    cmap.fun('color','Spectral',discrete=TRUE) +
    labs(y='Cases',x='Date') +
    theme_light()
  return(g)
}
plot.distr = function(C,y,ylabel,cmap='',...){
  g = ggplot(C,aes_string(y=y,...)) +
    geom_boxplot(alpha=.2,aes(x=decile)) +
    cmap.fun('color',cmap=cmap,discrete=TRUE) +
    cmap.fun('fill',cmap=cmap,discrete=TRUE) +
    labs(y=ylabel,x='Decile') +
    lims(y=c(0,30)) +
    theme_light()
  return(g)
}