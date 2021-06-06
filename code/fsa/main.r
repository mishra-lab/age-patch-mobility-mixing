args = commandArgs(trailingOnly=TRUE)

if (args[1] == 'epsilon'){
  source('epsilon.r')
  main.epsilon()
}
if (args[1] == 'mobility'){
  source('mobility.r')
  main.mobility()
}
if (args[1] == 'mixing'){
  source('mixing.r')
  for (mo in c('ref',mo.covid)){
    main.mixing(t=mo)
  }
  merge.save.mixing('Ci')
}
if (args[1] == 'debug'){
  # DEBUG
  source('runtime.r')
  # gen.mix.runtime.data()
  gen.Ci.gg.y()
}