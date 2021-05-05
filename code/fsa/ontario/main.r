source('data.r')
source('plot.r')

# b.s = 0.05 # >= max(rowSums(B.gg)) = 0.04017756s
b.h = 0.1 # proportion who mix at home with local traveller pool
b.y = c(
  'home'       = 0,
  'work'       = 1,
  'school'     = 0,
  'transport'  = 1,
  'leisure'    = .5,
  'otherplace' = .5)

B.ggo = load.fsa.mob()
B.gg = (1-b.h) * B.ggo / rowSums(B.ggo)
B.g  = rowSums(B.gg)
pop  = load.fsa.data()
pmc  = load.polymod.data()
C.ay = do.call(cbind,lapply(pmc,rowSums))
P.ga = P.n.to.g(pop)
Q.gy = P.ga %*% C.ay
G = length(info$decile)

B.gg.  = B.gg + diag(1-B.g)
T.gy   = t(B.gg.) %*% Q.gy
X.gg.y = list()
for (y in info$c.type){
  Q.gy. = B.gg. * Q.gy[,y]
  X.gg.y[[y]] = matrix(0,nrow=G,ncol=G)
  for (g in info$decile){
    X.gg.y[[y]] = X.gg.y[[y]] + outer(Q.gy.[,g],Q.gy.[,g]) / T.gy[g,y]
  }
  X.gg.y[[y]] = b.y[[y]] * X.gg.y[[y]] + diag((1-b.y[[y]]) * Q.gy[,y])
}

X.gg = Reduce('+',X.gg.y)
X.gg.y. = lapply(X.gg.y,function(x){x/1e6})
print(sapply(seq(6),function(y){ Q.gy[,y] / rowSums(X.gg.y[[y]]) })) # sum = 1
print(sapply(X.gg.y,function(x){ median(diag(x) / rowSums(x)) })) # % with self
print(median(diag(X.gg) / rowSums(X.gg)))

figname = function(f){ root.path('out','fig','fsa',paste0(f,'.pdf')) }

g = plot.mix(mix.melt(B.ggo),xy='decile')
ggsave(figname('mobility'),width=5,height=5)

g = plot.mix(mix.melt(X.gg/1e6),xy='decile')
ggsave(figname('Xgg'),width=5,height=5)

g = plot.mix(mix.melt(xfun(X.gg/1e6)),xy='decile')
ggsave(figname('Xgg-o'),width=5,height=5)

g = plot.mix.vs(X.gg.y.,clim=NULL,xy='decile')
ggsave(figname('Xggy'),width=14,height=3)

g = plot.mix.vs(xfun(X.gg.y.),clim=NULL,xy='decile')
ggsave(figname('Xggy-o'),width=14,height=3)

g = plot.mix.vs(pmc,clim=NULL,xy='age')
ggsave(figname('polymod'),width=14,height=3)