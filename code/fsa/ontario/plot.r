library(reshape2)
library(ggplot2)

.offd = function(x){ x-diag(rep(NA,nrow(x))) }
.norm = function(x){ x / max(x) }
.log  = function(x){ log(x + 1e-6*mean(x)) }

xfun = function(X,fun=.offd){
  if (is.list(X)){
    return (lapply(X,fun))
  } else {
    return (fun(X))
  }
}
mix.melt = function(X){
  return(melt(X,value.name='X',varnames=c('i','i.')))
}
plot.mix = function(X.,xy='age',clim=NULL,cmap='inferno'){
  g = ggplot(X.,aes(x=factor(i.),y=factor(i),fill=X,color=X)) +
    geom_tile() +
    coord_fixed(ratio=1) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    xlab(labels[[xy]]$x) + ylab(labels[[xy]]$y) +
    viridis::scale_fill_viridis(option=cmap,limits=clim) +
    viridis::scale_color_viridis(option=cmap,limits=clim) +
    theme_light() +
    guides(color=FALSE) +
    switch(xy,
      age = theme(axis.text.x=element_text(angle=90,hjust=1)),
      fsa = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    )
  return(g)
}
plot.mix.vs = function(X,xy='age',clim=NULL){
  ci.melt = function(group){
    ci = melt(X[[group]],value.name='X',varnames=c('i','i.'))
    ci$group = as.factor(group)
    return(ci)
  }
  X. = do.call(rbind,lapply(names(X),ci.melt))
  g = plot.mix(X.,xy=xy,clim=clim) +
    facet_grid(cols=vars(group)) +
    guides(fill=guide_colorbar(barheight=5))
  return(g)
}