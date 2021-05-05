source('data.r')
source('mix.r')
source('plot.r')

h.y = c(
  'home'       = 1,
  'work'       = 0,
  'school'     = 1,
  'transport'  = 0,
  'leisure'    = .5,
  'otherplace' = .5
)
pmc  = load.polymod.data()
pop  = load.fsa.data()
mob  = load.fsa.mob()
C.ay = pmc.to.Cay(pmc)
P.ga = pop.to.Pga(pop)
B.gg = mob$B.gg
X.gaga.y = gaga.mix(C.ay,P.ga,B.gg,h.y)

X.gaga = Reduce('+',X.gaga.y)
print(pct.self(a.sum(X.gaga,c(2,4))))
print(sapply(X.gaga.y,function(X){ pct.self(a.sum(X,c(2,4))) }))

figname = function(f){ root.path('out','fig','fsa',paste0(f,'.pdf')) }
g = plot.mix(B.gg,'g',aggr=FALSE);             ggsave(figname('mobility'),width= 5,height= 4)
g = plot.mix(X.gaga,'g');                      ggsave(figname('Xgg'),     width= 5,height= 4)
g = plot.mix(X.gaga,'g',xfun=offd);            ggsave(figname('Xgg-o'),   width= 5,height= 4)
g = plot.mix(X.gaga,'a');                      ggsave(figname('Xaa'),     width= 5,height= 4)
g = plot.mix(X.gaga,'a',xfun=offd);            ggsave(figname('Xaa-o'),   width= 5,height= 4)
g = plot.mix(X.gaga.y,'g');                    ggsave(figname('Xggy'),    width=14,height= 3)
g = plot.mix(X.gaga.y,'g',xfun=offd);          ggsave(figname('Xggy-o'),  width=14,height= 3)
g = plot.mix(X.gaga.y,'a',clim=c(0,10));       ggsave(figname('Xaay'),    width=14,height= 3)
g = plot.mix(pmc,'a',aggr=FALSE,clim=c(0,10)); ggsave(figname('polymod'), width=14,height= 3)
