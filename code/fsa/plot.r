suppressPackageStartupMessages({
library(reshape2)
library(viridis)
library(ggplot2)
library(ggridges)
})

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
plot.mix = function(C,what,vs,aggr=TRUE,xfun=NULL,clim=c(0,NA),cmap='inferno',...){
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
    scale_fill_viridis(option=cmap,limits=clim,end=.95,na.value='transparent') +
    scale_color_viridis(option=cmap,limits=clim,end=.95,na.value='transparent') +
    theme_light() +
    guides(color=FALSE,fill=guide_colorbar(barheight=5)) +
    switch(vs,
      a = theme(axis.text.x=element_text(angle=90,hjust=1)),
      n = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    )
  if (is.list(C)){
    if (length(C) <= 6) { g = g + facet_grid(cols=vars(group)) }
    else { g = g + facet_wrap(vars(group)) }
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
    g = g + guides(fill=FALSE)
  }
  return(g)
}