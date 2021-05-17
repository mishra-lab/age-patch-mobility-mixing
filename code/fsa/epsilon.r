library('reshape2')
library('ggplot2')
source('config.r')
source('plot.r')
source('mixing.r')
source('data.r')

# Objective: estimate epsilon for approximate representation of POLYMOD data on age mixing

plot.pmc.margins = function(pmc){
  margins = list(
    'Index' = rowSums,
    'Other' = colSums,
    'Ratio (Index / Other)' = function(x){ log(rowSums(x)/colSums(x)) }
  )
  c. = do.call(rbind,lapply(names(margins),function(m){
    c.m = sapply(pmc,margins[[m]])
    c.m. = mix.melt(c.m)
    c.m.$margin = m
    return(c.m.)
  }))
  # compare index / other margins: plot both, then ratio
  idx.mr = c.$margin == names(margins)[3]
  g = plot.mix(c.[!idx.mr,]) +
      facet_grid(cols=vars(margin)) + ylab('Age') + xlab('Contact Type')
  ggsave('../../out/fig/polymod-m.pdf',width=6,height=5)
  g = plot.mix(c.[idx.mr,],clim=c(-3.5,+3.5)) +
      facet_grid(cols=vars(margin)) + ylab('Age') + xlab('Contact Type')
  ggsave('../../out/fig/polymod-mr.pdf',width=4,height=5)
}
get.pmc.eps = function(pmc.y,eps.y){
  c.mi = rowSums(pmc.y)
  c.mo = colSums(pmc.y)
  c.r = c.mi %*% t(c.mo) / sum(c.mo)
  c.a = diag(c.mi)
  c.x = (1-eps.y) * c.r + (eps.y) * c.a
  rownames(c.x) = colnames(c.x)
  return(c.x)
}
get.pmc.error = function(pmc,pmc.eps){
  pmc = pmc + 1e-6
  return(sum(pmc * log(pmc / pmc.eps), na.rm=TRUE)) # KL-Divergence
  # return(mean(abs(pmc-pmc.eps)))     # absolute difference # DEBUG
  # return(mean(abs(pmc-pmc.eps)/pmc)) # relative difference # DEBUG
}
estimate.pmc.eps = function(pmc){
  eps.0 = c('home'=.5,'work'=.5,'school'=.5,'transport'=.5,'leisure'=.5,'otherplace'=.5)
  obj.fun = function(eps,aggr=TRUE) {
    pmc.eps = mapply(get.pmc.eps,pmc,eps,SIMPLIFY=FALSE)
    error = mapply(get.pmc.error,pmc,pmc.eps)
    if (aggr){ return(mean(error)) } else { return(error) }
  }
  opt = optim(eps.0,obj.fun,method='L-BFGS-B',lower=rep(0,6),upper=rep(1,6))
  eps.opt = opt$par
  print(eps.opt)
  print(obj.fun(eps.opt,aggr=FALSE))
  return(eps.opt)
}

main.epsilon = function(){
  pmc = load.polymod()
  eps = estimate.pmc.eps(pmc)
  pmc.eps = mapply(get.pmc.eps,pmc,eps,SIMPLIFY=FALSE)
  pmc.sum = lapply(list(
    'Original' = pmc,
    'Approx'   = pmc.eps
  ),function(x){ Reduce('+',x) })
  pmc.ratio = pmc.sum[[1]] / pmc.sum[[2]]

  g = plot.mix(pmc,    aggr=FALSE,clim=c(0, 5)); ggsave(figname('polymod','age-eps'),    width=14,height=3)
  g = plot.mix(pmc.eps,aggr=FALSE,clim=c(0, 5)); ggsave(figname('polymod-eps','age-eps'),width=14,height=3)
  g = plot.mix(pmc.sum,aggr=FALSE,clim=c(0,12)); ggsave(figname('polymod-vs','age-eps'), width= 5,height=3)
  g = plot.mix(pmc.ratio,aggr=FALSE,xfun=log10,clim=c(-1,+1),cmap='cividis');
    ggsave(figname('polymod-rs','age-eps'),width=4,height=3)
}






