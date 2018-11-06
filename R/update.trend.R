##############################################################
#' @title Recalculate posterior predictive trend with a different time scale.
#' @description Given a fitted model object (fit with 'keep.site.abund=TRUE') the function recalculates the posterior sample for an alternate
#' time scale than was originally specified in the model fitting MCMC sample.
#'
#' @param x An mcmc augmentation object produced by a call to \code{\link{mcmc.aggregate}} or
#' an element of the list produced by a call to \code{\link{newAggregation}}.
#' @param start A new start value for the time span
#' @param end A new end value for the time span
#' @param type The type of trend calculated. Use \code{"pred"} for posterior predictive trends
#' and \code{"real"} to use the estimated, realized abumndance aggregation.
#' @param order The order of trend calculated. Can be one of \code{"lin"}, for linear trends,
#' or, \code{"const"}, for mean log-abundence.
#'
#' @export
#'
updateTrend <- function(x, start, end, type="pred", order="lin"){
  #require(coda)
  if(type=="pred" | is.null(x$mcmc.sample$aggregated.real.abund)) smp <- x$mcmc.sample$aggregated.pred.abund
  else if(type=="real") smp <- x$mcmc.sample$aggregated.real.abund
  else stop("Unknown 'type', must be 'pred' or 'real'.\n")
  nms <- strsplit(colnames(smp),"-")
  time <- as.numeric(sapply(nms, function(x)x[[1]]))
  start <- max(start, min(time))
  end <- min(end, max(time))
  time.idx <- time>=start & time<=end
  agg <- sapply(nms, function(x)x[[2]])
  if(any(is.na(suppressWarnings(as.numeric(agg))))) agg <- factor(agg)
  else agg <- factor(as.numeric(agg))
  nms.agg <- as.character(levels(agg))

  yrs <- c(start:end)
  ag.df <- expand.grid(yrs, levels(agg))
  if(order=="lin"){
    if(length(nms.agg)>1) {
      ag.mm <- model.matrix(~(ag.df[,2]+0) + (ag.df[,1]:ag.df[,2]+0))
    } else {
      ag.mm <- cbind(rep(1,nrow(ag.df)), yrs)
    }
    H <- solve(crossprod(ag.mm))%*%t(ag.mm)
    y <- log(smp[,time.idx] + attr(x, "ln.adj"))
    tsmp <- t(apply(y, 1, function(y){H%*%y}))
    colnames(tsmp) <- c(paste(nms.agg, "(Intercept)"), paste(nms.agg, "(Trend)"))
  } else if(order=="const"){
    if(length(nms.agg)>1) ag.mm <- model.matrix(~(ag.df[,2]+0))
    else ag.mm <- matrix(rep(1,nrow(ag.df)), ncol=1)
    H <- solve(crossprod(ag.mm))%*%t(ag.mm)
    y <- log(smp[,time.idx] + attr(x, "ln.adj"))
    tsmp <- matrix(t(apply(y, 1, function(y){H%*%y})), ncol=ncol(ag.mm))
    colnames(tsmp) <- c(paste(nms.agg, "(Intercept)"))
  }
  else stop("Unknown 'type' specified! See ?updateTrend.")
  return(mcmc(tsmp))
}

