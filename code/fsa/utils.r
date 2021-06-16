# supporting functions for very general stuff

col.rename = function(x,old.name,new.name){
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
