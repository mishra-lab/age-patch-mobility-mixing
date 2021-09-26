suppressPackageStartupMessages({
  library(edfun)
})
source('utils.r')
source('data.r')
source('plot.r')

gen.YR = function(cases,rmin=1,rmax=3){
  Y = dcast(cases,t~decile)
  YR = cbind(Y[,2:10]/Y[,3:11],'10'=1)
  YR[YR<rmin] = rmin
  YR[YR>rmax] = rmax
  return(cbind.data.frame(YR,'t'=Y$t))
}

gen.YR.sampler = function(YR,trans=identity){
  return(apply(YR[,paste(seq(10))],2,function(x){
    suppressWarnings({
      if (any(x>1)){
        return( edfun(trans(x))$rfun )
      } else {
        return( function(n=1){ rep(1,n) } ) # TODO: trans?
      }
    })
  }))
}

apply.sampler = function(funs,n=1){
  return(sapply(funs,do.call,list(n=n)))
}

cum.YR = function(YR,...){
  cols = paste(seq(10))
  YRC = t(apply(YR[,rev(cols)],1,cumprod))[,cols]
  return(cbind.data.frame(YRC,...))
}

main.contacts = function(){
  cases = load.cases()
  YR = gen.YR(cases)
  trans = function(x){ return(c(x,runif(length(x),1,x))) }
  funs = gen.YR.sampler(YR,trans=trans)
  YRS  = apply.sampler(funs,nrow(YR))
  YRc  = cum.YR(YR,vs='Cases')
  YRSc = cum.YR(YRS,vs='Sample')
  ylabel = 'Ratio of Contacts vs Lowest Decile'
  Y.vs = melt(rbind(YRc,YRSc),variable.name='decile',value.name='ratio')
  YRC. = melt(YRc,variable.name='decile',value.name='ratio')
  plot.distr(Y.vs,'ratio',ylabel,color='vs',fill='vs')
    ggsave(figname('RCg-dvs','contacts'),width=6,height=4)
  plot.distr(YRC.,'ratio',ylabel,color='decile',fill='decile') +
    scale_deciles('color') + scale_deciles('fill') + theme(legend.position='none');
    ggsave(figname('RCg','contacts'),width=6,height=4)
  # print(colMeans(YRSc[,1:10])) # DEBUG
}