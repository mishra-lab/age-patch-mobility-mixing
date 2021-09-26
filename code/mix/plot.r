suppressPackageStartupMessages({
library(reshape2)
library(viridis)
library(ggplot2)
library(ggridges)
library(pals)
library(scales)
})
source('utils.r')

# supporting functions for plotting stuff

offd  = function(x){ x-diag(rep(NA,nrow(x))) }
list.bind = function(X,gname='group',fun=identity,arg.name=FALSE,...){
  # if X is a list, rbind its elements & add a new column bearing the list name
  # can run fun on each list element if needed, and possibly pass column name too
  if (is.list.(X)){
    return(do.call(rbind,lapply(names(X),function(name){
      if (arg.name){
        Xi = fun(X[[name]],name,...)
      } else {
        Xi = fun(X[[name]])
      }
      Xi[[gname]] = factor(name,levels=names(X))
      return(Xi)
    })))
  } else {
    return(fun(X))
  }
}
mix.melt = function(C,what,vs,aggr=TRUE,xfun=NULL,...){
  # list.bind + need to aggregate & melt each individual matrix
  if (is.null(xfun)){ xfun = identity }
  return(list.bind(C,gname='group',fun=function(Ci){
    return(melt(xfun(aggr.mix(Ci,what,vs,aggr=aggr,...)),
      value.name='X',varnames=c('i','i.')))
  }))
}
nsqrt_trans = function(){ # sqrt transform that preserves -/+
  trans_new('nsqrt',function(x){sign(x)*sqrt(abs(x))},function(x){sign(x)*x^2})
}
scale_deciles = function(what='color',...){
  return(switch(what,
    color = scale_color_discrete(type=rev(tol.rainbow(10)),...),
    fill  = scale_fill_discrete( type=rev(tol.rainbow(10)),...)
  ))
}
plot.mix = function(C,what,vs,aggr=FALSE,xfun=NULL,
    trans='sqrt',clim=c(NA,NA),cmap='inferno',gez=TRUE,...){
  C. = mix.melt(C,what,vs,aggr=aggr,xfun=xfun,...)
  if (gez){ C.$X[C.$X<0] = 0 }
  v.name = switch(what,
    CX = 'Total contacts\n(Millions)',
    Ci = 'Contacts\nper person',
    Cp = 'Contact\nformation\nprobability',
    B  = '% decile g\nwho are mobile\nin decile g\'',
    Bc = '% decile g\ntravellers who\nvisit decile g\'')
  g = ggplot(C.,aes(x=factor(i),y=factor(i.),fill=X,color=X)) +
    geom_tile() +
    coord_fixed(ratio=1) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    labs(x=config$labels[[vs]]$x,y=config$labels[[vs]]$y,fill=v.name) +
    theme_light() +
    guides(color='none',fill=guide_colorbar(barheight=5)) +
    switch(vs,
      a = theme(axis.text.x=element_text(angle=90,hjust=1)),
      n = theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    ) +
    switch(cmap,
      inferno = scale_fill_viridis_c(option='B',end=.95,limits=clim,trans=trans),
      RdBu    = scale_fill_distiller(palette='RdBu',limits=clim,trans=trans)) +
    switch(cmap,
      inferno = scale_color_viridis_c(option='B',end=.95,limits=clim,trans=trans),
      RdBu    = scale_color_distiller(palette='RdBu',limits=clim,trans=trans))
  if (is.list.(C)){
    if (length(C) <= 6) { # HACK
      g = g + facet_grid(cols=vars(group))
    } else {
      g = g + facet_wrap(vars(group),ncol=6) +
        theme(legend.position=c(.84,-.03),legend.justification=c(0,0))
    }
  }
  return(g)
}
void = function(){
  return(theme_void() +
    theme(strip.text.x=element_blank(),strip.text.y=element_blank()))
}
plot.ridge.density = function(X,x,y='month',fill='q4',xmax=NULL,legend=FALSE){
  if (is.null(xmax)){ xmax=max(X[[x]],na.rm=TRUE) }
  if (y=='month'){ y.lab='Month' } else { y.lab = y }
  X. = list.bind(X)
  g = ggplot(X.,aes_string(y=y,x=x))
  if (fill=='q4'){
    g = g + stat_density_ridges(
        aes_string(fill='factor(stat(quantile))'),
        geom='density_ridges_gradient',calc_ecdf=TRUE,quantiles=4) +
      scale_fill_viridis_d(option='B',alpha=.7,begin=.2,end=.8)
    }
  if (fill=='decile'){
    g = g + stat_density_ridges(aes_string(fill='decile'),alpha=.5,color=NA) +
      scale_deciles('fill')
  }
  g = g + scale_y_discrete(limits=rev) +
    xlim(0,xmax) + labs(x=x,y=y.lab,fill=fill) +
    theme_light()
  if (!legend){
    g = g + guides(fill='none')
  }
  if (is.list.(X)) { g = g + facet_grid(cols=vars(group)) }
  return(g)
}
plot.pop.density = function(){
  pop = load.fsa.pop(aggr=FALSE)
  P = aggregate(pop~group+age,pop,sum)
  g = ggplot(P,aes(x=age,y=pop/mean(pop))) +
    geom_line(aes(color=group)) +
    scale_deciles('color') +
    labs(x='Age',y='Relative Population Proportion',color='Decile') +
    theme_light()
  return(g)
}
plot.contact.margins = function(C.aa.y.list,age.list,style='line'){
  C.a.y. = list.bind(C.aa.y.list,gname='vs',arg.name=TRUE,function(C.aa.y,C.name){
    list.bind(C.aa.y,gname='Type',function(C.aa){
      return(data.frame(Contacts=rowSums(C.aa),Age=midpoint(age.list[[C.name]])))
    })
  })
  g = ggplot(C.a.y.,aes(x=Age,y=Contacts,group=vs,color=vs)) +
    facet_grid(cols=vars(Type)) +
    scale_color_viridis_d(option='B',begin=.15,end=.85) +
    theme_light()
  if (style=='line'){ g = g + geom_line() }
  if (style=='dots'){ g = g + geom_line() + geom_point(size=1,alpha=.5) }
  if (style=='area'){ g = g + geom_area(aes(fill=vs),position='identity',alpha=.5) }
  if (style=='bars'){ g = g + geom_bar(aes(fill=vs),stat='identity',position='identity',alpha=.5) }
  return(g)
}
plot.cond.mob = function(Xgg,style='line'){
  g = ggplot(Xgg,aes(y=Bc,color=decile.visited,group=decile.visited)) +
    facet_wrap(vars(decile),ncol=5) +
    scale_deciles('color') +
    labs(y='Conditional destination probability') +
    theme_light() +
    guides(color=guide_legend(ncol=10)) + theme(legend.position='top')
  if (style=='line'){
    g = g + geom_line(aes(x=month)) + xlab('Month') +
      scale_x_discrete(guide = guide_axis(angle = 90)) }
  if (style=='box') {
    g = g + geom_boxplot(aes(x=decile.visited)) + xlab('Decile Visited') }
  return(g)
}
plot.p.mob = function(Tg,y,style='line'){
  y.lab = switch(y,
    'rho.gt/rho.t' = 'Relative time away from home vs mean',
    'rho.t'        = 'Mean time away from home',
    'phi.gt/phi.t' = 'Relative proportion of away time within home FSA vs mean',
    'phi.t'        = 'Mean proportion of away time within home FSA',
  )
  g = ggplot(Tg,aes_string(y=y,color='decile',group='decile')) +
    scale_deciles('color') +
    ylab(y.lab) +
    theme_light()
  if (style=='line'){
    g = g + geom_line(aes(x=ordered(month))) + xlab('Month') +
      scale_x_discrete(guide = guide_axis(angle = 90)) }
  if (style=='mean'){
    g = g + geom_line(aes(color=NULL,x=ordered(month))) + xlab('Month') +
      scale_x_discrete(guide = guide_axis(angle = 90)) }
  if (style=='box'){
    g = g + geom_boxplot(aes(x=decile)) + xlab('Home Decile') }
  return(g)
}
plot.B.approx = function(Xgg){
  g = ggplot(Xgg,aes(x=B,y=Ba,color=decile.visited)) +
    geom_point(size=0.5) +
    facet_wrap(vars(decile),ncol=5) +
    scale_x_continuous(trans='log10') + xlab('Observed Mobility') +
    scale_y_continuous(trans='log10') + ylab('Approximated Mobility') +
    scale_deciles('color') +
    theme_light() +
    guides(color=guide_legend(ncol=10)) + theme(legend.position='top')
  return(g)
}
plot.cases = function(cases){
  g = ggplot(cases,aes(x=as.Date(t),y=cases,group=decile,color=decile)) +
    geom_line() +
    scale_deciles('color') +
    labs(y='Cases',x='Date') +
    theme_light()
  return(g)
}
plot.distr = function(C,y,ylabel,...){
  g = ggplot(C,aes_string(y=y,...)) +
    geom_boxplot(alpha=.2,aes(x=decile)) +
    labs(y=ylabel,x='Decile') +
    lims(y=c(0,30)) +
    theme_light()
  return(g)
}
plot.map = function(shp=NULL){
  if (is.null(shp)){ shp = load.shape() }
  shp$decile = factor(shp$decile,levels=seq(10))
  g = ggplot(shp) +
    geom_sf(aes(fill=decile),lwd=.1,color='white') +
    scale_deciles('fill',drop=FALSE) +
    labs(fill='Decile') +
    theme_light()
  return(g)
}
plot.map.main = function(ext='.png'){
  shp = load.shape()
  shp = merge(shp,read.csv(root.path('data','mix','fsa_region.csv')))
  plot.map(subset(shp,substr(FSA,1,1)=='P')) + theme_void() + guides(fill='none');
    ggsave(figname('ontario-north','map',ext=ext),w=3,h=3,dpi=600)
  plot.map(subset(shp,substr(FSA,1,1)!='P')) + theme_void() + guides(fill='none');
    ggsave(figname('ontario-south','map',ext=ext),w=3,h=3,dpi=600)
  plot.map(subset(shp,GTA==1)) + theme_void();
    ggsave(figname('ontario-gta','map',ext=ext),w=4,h=3,dpi=600)
}
