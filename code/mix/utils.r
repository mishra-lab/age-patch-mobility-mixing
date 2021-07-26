# supporting functions for very general stuff

col.rename = function(x,old.name,new.name){ # rename a column
  colnames(x)[grepl(old.name,colnames(x))] = new.name
  return(x)
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
midpoint = function(x,ends=FALSE){
  # get the midpoints of vector x, and repeat the last value
  return(x+diff(bin.extrap(x))/2)
}
bin.extrap = function(x,post=1,pre=0,x.max=NULL,x.min=NULL){
  # e.g. if x = [0,5,10,15] represents 0-4, 5-9, 10-14, 15-19.
  # then we extrapolate additional bins either n bins = pre/post,
  # or up to x.min/x.max, however many are needed
  l = length(x)
  d.0 = x[2] - x[1]
  d.1 = x[l] - x[l-1]
  if (!is.null(x.min)) { x.pre  = seq(x[1]-d.0, x.min, -d.0) }
  else if (pre > 0)    { x.pre  = seq(x[1]-d.0, x[1]-d.0*pre, -d.0) }
  else                 { x.pre  = numeric(0) }
  if (!is.null(x.max)) { x.post = seq(x[l]+d.1, x.max, +d.1) }
  else if (post > 0)   { x.post = seq(x[l]+d.1, x[l]+d.1*post, +d.1) }
  else                 { x.post = numeric(0) }
  return(c(rev(x.pre),x,x.post))
}
bin.size = function(x,norm=FALSE){
  # get widths of a stratification "x", and possibly normalize.
  # e.g. if [x] is [1,3,7], bs is [2,4,4], and if norm, bs is [.2,.4,.4]
  # since bin.extrap repeats the last width (7-3=4)
  bs = diff(bin.extrap(x))
  names(bs) = names(x)
  if (norm){ bs = bs/sum(bs) }
  return(bs)
}