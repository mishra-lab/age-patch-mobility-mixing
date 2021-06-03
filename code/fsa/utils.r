col.rename = function(x,old.name,new.name){
  colnames(x)[grepl(old.name,colnames(x))] = new.name
  return(x)
}
a.sum = function(A,d){
  ds = seq(length(dim(A)))
  return(colSums(aperm(A,c(d,ds[-d])),dims=length(d)))
}
a.mean = function(A,d,W=NULL){
  if (is.null(W)){ W = array(1,dim(A)) }
  ds = seq(length(dim(A)))
  dp = c(d,ds[-d])
  return(colSums(aperm(A*W,dp),dims=length(d)) / colSums(aperm(W,dp),dims=length(d)))
}
interp1d = function(yi,xi,xo){
  return(approx(xi,yi,xo,method='linear',yleft=yi[1],yright=yi[length(yi)])$y)
}
interp2d = function(Y,xi,xo){
  return(apply(apply(Y,1,interp1d,xi=xi,xo=xo),1,interp1d,xi=xi,xo=xo))
}