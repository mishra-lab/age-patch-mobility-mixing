source('../config.r')

load.fsa.data = function(){
  pop = read.csv(root.path('data','fsa','population.csv'))
  return(pop)
}

load.fsa.mob = function(){
  mob = read.csv(root.path('data','fsa','mobility.csv'))
  B = matrix(mob$visiting_mean,nrow=10,ncol=10)
  colnames(B) = info$decile
  rownames(B) = info$decile
  return(list(B.gg=B))
}

load.polymod.data = function(){
  pmc = list()
  f = 1
  for (c.type in info$c.type){
    names(f) = paste('cnt',c.type,sep='_')
    pmc.y = suppressMessages(contact_matrix(
      polymod,filter=f,
      age.limits=info$age
      # age.limits=c(0,20,40,60,80) # DEBUG
      # age.limits=c(0,15,30,45,60,75) # DEBUG
      # age.limits=c(0,10,20,30,40,50,60,70,80) # DEBUG
    ))$matrix
    rownames(pmc.y) = names(info$age)
    colnames(pmc.y) = names(info$age)
    pmc[[c.type]] = pmc.y
  }
  return(pmc)
}
