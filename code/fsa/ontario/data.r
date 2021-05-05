source('../config.r')

load.fsa.data = function(){
  pop = read.csv(root.path('data','population.csv'))
  return(pop)
}

load.fsa.mob = function(){
  mob = read.csv(root.path('data','mobility.csv'))
  B = matrix(mob$visiting_mean,nrow=10,ncol=10)
  colnames(B) = info$decile
  rownames(B) = info$decile
  return(B)
}

load.polymod.data = function(){
  pmc = list()
  f = 1
  for (c.type in info$c.type){
    names(f) = paste('cnt',c.type,sep='_')
    pmc[[c.type]] = suppressMessages(contact_matrix(
      polymod,filter=f,
      age.limits=info$age
      # age.limits=c(0,20,40,60,80) # DEBUG
      # age.limits=c(0,15,30,45,60,75) # DEBUG
      # age.limits=c(0,10,20,30,40,50,60,70,80) # DEBUG
    ))$matrix
  }
  return(pmc)
}

P.n.to.g = function(pop){
  G.n  = matrix(pop$decile,nrow=length(info$age))[1,] # WARNING: not robust
  S.gn = outer(info$decile,G.n,'==')*1
  P.na = matrix(pop$pop_agefsa,nrow=length(info$age)) # WARNING: not robust
  P.ga = S.gn %*% t(P.na)
  colnames(P.ga) = names(info$age)
  return(P.ga)
}