get_mode = function(x){
  ny <- length(x)
  k <- ceiling(ny/2) - 1
  y <- sort(x)
  inf <- y[1:(ny-k)]
  sup <- y[(k+1):ny]
  diffs <- sup - inf
  i <- which(diffs==min(diffs))
  return(mean(y[i:(i+k)]))
}
