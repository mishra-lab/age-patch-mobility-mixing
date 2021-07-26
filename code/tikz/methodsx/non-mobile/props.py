import numpy as np
_ = None

# compute M_ij: the proportions of i's contacts formed among j
# under 2 different interpretations of non-mobile:
# "home": non-mobile don't mix with other travellers
# "stay": non-mobile stay in their patch but may mix with other travellers
# optionally, weight the "home" & "travel" contacts by ch & ct, resp
# -- note: in case=="stay", ch applies to non-mobile,
#          even though they mix with mobile (as in Arenas 2020)

def mix(mi,case,ch=1,ct=1,nmh=False):
  n = len(mi)
  pi = np.ones(n)
  mi = np.array(mi,dtype=float)
  Bii = np.ones((n,n))/n
  Mii = np.zeros((n,n))
  for x in range(n): # for each pool
    Qx = pi * mi * Bii[:,x] * (ct + (ch if nmh else 0)) # mobile pop in this pool
    if case == 'home':
      Mii[x,x] += pi[x] * (1-mi[x]) * (ch) # non-mobile pop: add 100% on diag
    if case == 'stay':
      Qx[x] += pi[x] * (1-mi[x]) * (ch) # non-mobile pop: add to this pool
    X = np.outer(Qx,Qx/Qx.sum()) # proportionate mixing in this pool
    X[np.isnan(X)] = 0
    Mii += X # add across pools
  return Mii/Mii.sum(axis=1,keepdims=True) # we want % contacts, not abs counts

def pext(Mii):
  return (1-np.diag(Mii))

print(pext(mix([1. ,1. ,1. ],'home')))
print(pext(mix([1. ,1. ,1. ],'stay')))
print(pext(mix([ .0, .0, .0],'home')))
print(pext(mix([ .0, .0, .0],'stay')))
print(pext(mix([ .5, .5, .5],'home',nmh=True)))
print(pext(mix([ .5, .5, .5],'stay',nmh=True)))
print(pext(mix([ .0, .5,1. ],'home',nmh=True)))
print(pext(mix([ .0, .5,1. ],'stay',nmh=True)))

