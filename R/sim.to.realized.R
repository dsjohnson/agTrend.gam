
sim.to.realized = function(data, post.sample){
  out=post.sample
  for(j in 1:nrow(data)){
    out[data$time[j]+1,] = data$y[j]
  }
  return(out)
}
