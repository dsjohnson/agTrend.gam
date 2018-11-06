#######################################################################
#' @title Fit a gam model to an individual site data set
#' @description The purpose of this function is to allow a user to fit the gam model which is used for
#' imputation to a sigle site's worth of data. This is primarily used for testing purposes and not meant
#' to be called by the user in general.
#' @param data data
#' @param obs.formula Formula for the observation model
#' @param sig.abund Variance estimates for abundance response if it is estimated.
#' @param min.k Minimum df for the thin-plate regression spine model fit to the counts
#' @param ln.adj Amount to add to zero valued abundances before taking logs
#' @param obs.prior List that provides the mean (named 'gamma.mean') and covariance matrix
#' for the observation model coefficients (named 'gamma.vcov')
#' @return a list...
#' @author Devin S. Johnson
#' @export
#' @import mgcv
site.fit.gam = function(data, obs.formula=NULL, sig.abund=NULL,
                        min.k = 3, ln.adj=0, family="tw",
                        obs.prior=NULL){
  nz = sum(data$y>0)

  if(!is.null(obs.formula)){
    if(is.null(obs.prior)) stop("A list of priors must be given if 'obs.formula' is specified!")
    # gamma offset/random eff.
    X = model.matrix(obs.formula, data)
    Xgamma = X%*% obs.prior$gamma.mean %>% as.vector()
    gamma_pen = list(X=list(S=solve(obs.prior$gamma.vcov), sp=1))
    gam.cut = ncol(X)+min.k+3
    lin.cut = ncol(X)+3
    if(all(apply(X,2,var)==0)){
      use.gamma=FALSE
      } else{use.gamma=TRUE}
  } else(us.gamma=FALSE)
  if(use.gamma){
    # Fit models
    if(nz>=gam.cut){
      fit = mgcv::gam(
        y ~ offset(Xgamma) + X + s(time, k=nrow(data)-ncol(X)-1),
        select=TRUE, data=data, family=family, method="REML",
        paraPen = gamma_pen
      )
    } else if(nz>=lin.cut){
      pen_list = c(gamma_pen, list(time=list(sp=-1, S=diag(1))))
      fit <- mgcv::gam(
        y ~ offset(Xgamma) + time + X,
        paraPen = pen_list, data=data, family=family, method="REML"
      )
    } else{
      fit <- mgcv::gam(
        y ~ offset(Xgamma) + X,
        paraPen = gamma_pen, data=data, family=family, method="REML"
      )
    }
  } else{
    # Fit models
    if(nz>=min.k+3){
      fit = mgcv::gam(
        y ~ s(time, k=nrow(data)-1),
        select=TRUE, data=data, family=family, method="REML"
      )
    } else if(nz>= 3){
      pen_list = list(time=list(sp=-1, S=diag(1)))
      fit = mgcv::gam(
        y ~ time,
        paraPen = pen_list, data=data, family=family, method="REML"
      )
    }
    else {
      fit = mgcv::gam(
        y ~ 1,
        data=data, family=family, method="REML"
      )
    }
  } # end if obs.formula present
  return(fit)
}
