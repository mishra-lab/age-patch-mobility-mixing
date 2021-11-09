args = commandArgs(trailingOnly=TRUE)
source('config.r')
set.config()

if (args[1] == 'mobility'){
  source('mobility.r')
  main.mobility()
}
if (args[1] == 'mixing'){
  source('mixing.r')
  main.mixing()
}
if (args[1] == 'assumptions'){
  source('mixing.r')
  mixing.assumptions()
}
if (args[1] == 'debug'){
  source('mixing.r')
  # gen.save.mixing(n.y=4)
  plot.map.main()
  plot.cases(load.cases());
  ggsave(figname('Yg'),w=8,h=4)
}