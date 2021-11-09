# supporting functions for very general stuff

col.rename = function(x,old.name,new.name){ # rename a column
  colnames(x)[grepl(old.name,colnames(x))] = new.name
  return(x)
}
is.list. = function(x){
  return(class(x)=='list')
  # return(inherits(x,'list'))
}
a.sum = function(A,d){ # sum across dimensions "d" of N-D array "A"
  if (length(d)==0){ return(A) }
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}
a.mean = function(A,d,W=NULL){ # mean across dimensions "d" of N-D array "A"
  # possibly weighted by "W" having same dim as "A"
  if (length(d)==0){ return(A) }
  if (is.null(W)){ W = array(1,dim(A)) }
  ds = seq(length(dim(A)))
  dp = c(d,ds[-d])
  return(colSums(aperm(A*W,dp),dims=length(d)) / colSums(aperm(W,dp),dims=length(d)))
}
symmetric = function(X){
  # force X to be symmetric by averaging with the transpose
  # X can be 2D, 4D, 6D, etc. but the dimensions must be ordered i1,i2,...j1,j2,...
  n = length(dim(X))
  d = seq(n)
  return((X+aperm(X,c(d[(n/2+1):n],d[1:(n/2)])))/2)
}
interp2d = function(Y,xi,xo){
  xi[length(xi)] = xi[length(xi)] - 1e-9 # HACK akima is broken IRDFK
  g = expand.grid(x=xi,y=xi)
  Yo = akima::interp(x=g$x, y=g$y, z=as.vector(Y), xo=xo, yo=xo, linear=TRUE)$z
  return(Yo)
}
diag.pad = function(X,pre=0,post=0){
  # diagonally pad matrix X by pre (i < 0) & post (i > nrow(X)); assume X is square
  N0 = nrow(X)
  N1 = N0 + pre + post
  X1 = matrix(0,nrow=N1,ncol=N1)
  X1[pre+(1:N0),pre+(1:N0)] = X
  X1[1:pre,pre+N0+1:post] = X[1,N0]
  X1[pre+N0+1:post,1:pre] = X[N0,1]
  for (i in seq(pre,1)){
    X1[i,i-1+(1:N0)] = X[1,1:N0]
    X1[i,N0+(i:pre)] = X[1,N0]
    X1[i-1+(1:N0),i] = X[1:N0,1]
    X1[N0+(i:pre),i] = X[N0,1]
  }
  for (i in seq(1,post)){
    X1[pre+N0+i,pre+i+(1:N0)] = X[N0,1:N0]
    X1[pre+N0+i,pre+(1:i)]    = X[N0,1]
    X1[pre+i+(1:N0),pre+N0+i] = X[1:N0,N0]
    X1[pre+(1:i),pre+N0+i]    = X[1,N0]
  }
  return(X1)
}
lm1 = function(x){
  return(x[1:(length(x)-1)])
}
age.names = function(age){
  return(lm1(names(age)))
}
midpoint = function(x){
  # get the midpoints of vector x, and repeat the last value
  return(lm1(x)+diff(x)/2)
}
bin.size = function(x,norm=FALSE){
  # get widths of a stratification "x", and possibly normalize.
  # e.g. if x & x.max are [1,3,7] & 11, bs is [2,4,4], and if norm, bs is [.2,.4,.4]
  bs = diff(x)
  names(bs) = lm1(names(x))
  if (norm){ bs = bs/sum(bs) }
  return(bs)
}