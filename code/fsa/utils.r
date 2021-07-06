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
interp1d = function(yi,xi,xo){
  # interpolate yi from points xi to xo (linear with repeated endpoints)
  return(approx(xi,yi,xo,method='linear',yleft=yi[1],yright=yi[length(yi)])$y)
}
interp2d = function(Y,xi,xo){
  # interpolate Y (square) from points xi to xo in both dimensions
  return(apply(apply(Y,1,interp1d,xi=xi,xo=xo),1,interp1d,xi=xi,xo=xo))
}
midpoint = function(x){
  # get the midpoints of vector x, and repeat the last value
  return(x+c(diff(x)/2,0))
}
bin.size = function(x,x.max,norm=FALSE){
  # get widths of a stratification "x", and possibly normalize.
  # e.g. if [x,x.max] is [1,2,5,11], bs is [1,3,6], and if norm, bs is [.1,.3,.6]
  bs = diff(c(x,x.max))
  names(bs) = names(x)
  if (norm){ bs = bs/sum(bs) }
  return(bs)
}