library(reshape2)
library(ggplot2)

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
    viridis::scale_fill_viridis(option=cmap,limits=clim) +
    viridis::scale_color_viridis(option=cmap,limits=clim) +
    theme_light() +
    guides(color=FALSE,fill=guide_colorbar(barheight=5)) +
    switch(case,
      a = theme(axis.text.x=element_text(angle=90,hjust=1)),
      n = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    )
  if (is.list(X)){ g = g + facet_grid(cols=vars(group)) }
  return(g)
}