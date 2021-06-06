library('reshape2')
library('ggplot2')
source('config.r')
source('plot.r')
source('mixing.r')
source('data.r')

# Objective: estimate epsilon for approximate representation of POLYMOD data on age mixing

plot.cy.margins = function(C.y){
  margins = list(
    'Index' = rowSums,
    'Other' = colSums,
    'Ratio (Index / Other)' = function(x){ log(rowSums(x)/colSums(x)) }
  )
  c. = do.call(rbind,lapply(names(margins),function(m){
    c.m = sapply(C.y,margins[[m]])
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
get.cy.eps = function(C.yi,eps.y){
  c.mi = rowSums(C.yi)
  c.mo = colSums(C.yi)
  c.r = c.mi %*% t(c.mo) / sum(c.mo)
  c.a = diag(c.mi)
  c.x = (1-eps.y) * c.r + (eps.y) * c.a
  rownames(c.x) = colnames(c.x)
  return(c.x)
}
get.cy.error = function(C.yi,C.yi.eps){
  C.yi = C.yi + 1e-6
  return(sum(C.yi * log(C.yi / C.yi.eps), na.rm=TRUE)) # KL-Divergence
  # return(mean(abs(C.yi-C.yi.eps)))      # absolute difference # DEBUG
  # return(mean(abs(C.yi-C.yi.eps)/C.yi)) # relative difference # DEBUG
}
estimate.cy.eps = function(C.y){
  eps.0 = c('home'=.5,'other'=.5)
  obj.fun = function(eps,aggr=TRUE) {
    C.y.eps = mapply(get.cy.eps,C.y,eps,SIMPLIFY=FALSE)
    error = mapply(get.cy.error,C.y,C.y.eps)
    if (aggr){ return(mean(error)) } else { return(error) }
  }
  opt = optim(eps.0,obj.fun,method='L-BFGS-B',lower=rep(0,6),upper=rep(1,6))
  eps.opt = opt$par
  print(eps.opt)
  print(obj.fun(eps.opt,aggr=FALSE))
  return(eps.opt)
}

main.epsilon = function(){
  C.y = load.contacts()
  eps = estimate.cy.eps(C.y)
  C.y.eps = mapply(get.cy.eps,C.y,eps,SIMPLIFY=FALSE)
  C.y.sum = lapply(list(
    'Original' = C.y,
    'Approx'   = C.y.eps
  ),function(x){ Reduce('+',x) })
  C.y.ratio = C.y.sum[[1]] / C.y.sum[[2]]

  f = 'epsilon'
  plot.mix(C.y,    aggr=FALSE,clim=c(0,10)); ggsave(figname('canada',    f), width=8,height=4)
  plot.mix(C.y.eps,aggr=FALSE,clim=c(0,10)); ggsave(figname('canada-eps',f), width=8,height=4)
  plot.mix(C.y.sum,aggr=FALSE,clim=c(0,10)); ggsave(figname('canada-vs', f), width=8,height=4)
  plot.mix(C.y.ratio,aggr=FALSE,xfun=log10,clim=c(-1,+1),cmap='cividis'); ggsave(figname('canadaf-rs',f),width=5,height=4)
}

