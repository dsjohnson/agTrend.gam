#######################################################################
#' @title Simulate abundance measurements from a fitted GAM model
#' @description The purpose of this function is to allow a user to simulate abundance
#' measurements from a fitted GAM model. This is primarily used for testing purposes and not meant
#' to be called by the user in general.
#' @param fits Fitted gam object from mgcv
#' @param data data
#' @param max.time The maximum time used for simulation
#' @param size Size of the drawn sample
#' @return a matrix that is max.time by size
#' @author Devin S. Johnson
#' @export
#' @import mgcv tweedie
#' @import mvtnorm

site.sim.gam = function(fits, data, max.time, size){
  b = coef(fits)
  V = vcov(fits, unconditional = TRUE)
  newdata = list(
    time=0:max.time,
    Xgamma=rep(0, max.time+1),
    X = matrix(0, nrow=max.time+1, ncol=max(1,length(grep("X", names(b)))))
  )
  if(length(b)>1){
    L = mgcv::predict.gam(fits, newdata=newdata, type="lpmatrix")
  } else{
    L = matrix(1, nrow=length(newdata$time), ncol=1)
  }
  b.sim = mvtnorm::rmvnorm(size, b, V)
  mu.sim = L %*% t(b.sim) %>% base::exp() %>% as.vector()
  p = fits$family$getTheta(trans = T)
  phi = summary(fits)$dispersion
  a.sim = tweedie::rtweedie(n=length(mu.sim), mu=mu.sim, phi=phi, power=p) %>% matrix(nrow=nrow(L),ncol=size) %>% round()
  return(a.sim)
}
