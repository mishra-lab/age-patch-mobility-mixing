suppressPackageStartupMessages({
library(reshape2)
library(viridis)
library(ggplot2)
library(ggridges)
})

offd  = function(x){ x-diag(rep(NA,nrow(x))) }

mix.melt = function(X,aggr,xfun=NULL){
  if (isFALSE(aggr)){ aggr.fun = function(X){ X }}
  if (aggr=='a')    { aggr.fun = function(X){ a.sum(X/1e6,c(1,3)) }}
  if (aggr=='g')    { aggr.fun = function(X){ a.sum(X/1e6,c(2,4)) }}
  if (is.null(xfun)){ xfun = identity }
  if (is.list(X)){
    return(do.call(rbind,lapply(names(X),function(group){
      x = melt(xfun(aggr.fun(X[[group]])),value.name='X',varnames=c('i','i.'))
      x$group = factor(group,levels=names(X))
      return(x)
    })))
  } else {
    return( melt(xfun(aggr.fun(X)),value.name='X',varnames=c('i','i.')) )
  }
}
plot.mix = function(X,case='a',xfun=NULL,aggr=NULL,clim=NULL,cmap='inferno'){
  aggr = ifelse(is.null(aggr),case,aggr)
  X. = mix.melt(X,aggr,xfun)
  g = ggplot(X.,aes(x=factor(i.),y=factor(i),fill=X,color=X)) +
    geom_tile() +
    coord_fixed(ratio=1) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    xlab(labels[[case]]$x) + ylab(labels[[case]]$y) +
    scale_fill_viridis(option=cmap,limits=clim,na.value='transparent') +
    scale_color_viridis(option=cmap,limits=clim,na.value='transparent') +
    theme_light() +
    guides(color=FALSE,fill=guide_colorbar(barheight=5)) +
    switch(case,
      a = theme(axis.text.x=element_text(angle=90,hjust=1)),
      n = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    )
  if (is.list(X)){
    if (length(X) <= 6) { g = g + facet_grid(cols=vars(group)) }
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