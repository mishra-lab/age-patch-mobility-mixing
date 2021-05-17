args = commandArgs(trailingOnly=TRUE)

if (args[1] == 'epsilon'){
  source('epsilon.r')
  main.epsilon()
}
if (args[1] == 'mobility'){
  source('mobility.r')
  main.mobility()
}
if (args[1] == 'mixing')  {
  source('mixing.r')
  main.mixing()
}
