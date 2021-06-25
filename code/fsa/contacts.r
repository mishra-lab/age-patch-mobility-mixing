suppressPackageStartupMessages({
  library(edfun)
})
source('data.r')
source('plot.r')

gen.YR = function(cases,rmin=1,rmax=3){
  Y = dcast(cases,t~decile)
  YR = cbind(Y[,2:10]/Y[,3:11],'10'=1)
  YR[YR<rmin] = rmin
  YR[YR>rmax] = rmax
  return(cbind.data.frame(YR,'t'=Y$t))
}

gen.YR.sampler = function(YR){
  return(apply(YR[,paste(seq(10))],2,function(x){
    suppressWarnings({
      if (any(x>1)){
        return( edfun(x)$rfun )
      } else {
        return( function(n=1){ rep(1,n) } )
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
  funs = gen.YR.sampler(YR)
  YRS  = apply.sampler(funs,nrow(YR))
  YRc  = cum.YR(YR,vs='Data')
  YRSc = cum.YR(YRS,vs='Sample')
  y = melt(rbind(YRc,YRSc),variable.name='decile',value.name='ratio')
  ylabel = 'Ratio of Cases vs Lowest Decile'
  g = plot.distr(y,'ratio',ylabel,color='vs',fill='vs'); ggsave(figname('RCg','contacts'),width=6,height=4)
}