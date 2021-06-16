args = commandArgs(trailingOnly=TRUE)
source('config.r')
set.config('10x10')

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
  for (t in c('REF',config$t.covid)){
    main.mixing(t=t);
  }
  merge.save.mixing('Ci')
}
if (args[1] == 'debug'){
  # DEBUG
  source('runtime.r')
  test.runtime('2x2')
  test.runtime('10x10')
}
