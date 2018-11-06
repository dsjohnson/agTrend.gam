##############################################################
#' @title Recalculate new site aggregations from a previous MCMC aggregation
#'
#' @description
#' If the site abundence sample was retained in a call to \code{\link{mcmc.aggregate}} the
#' sites can be re-aggregated according to different region specifications.
#'
#' @param fit The output list from a previous call to \code{\link{mcmc.aggregate}}.
#' In order to use this function, \code{keep.site.abund = TRUE} had to be used in the original creation of \code{fit}. Else,
#' there is nothing to be aggregated!
#' @param aggregation.data  A data frame with the sites in one column
#' (with the same name as \code{site.name} used in the original call to create \code{fit}). The other columns
#' are factor variables defining other site aggregations.
#' @param type Which site abundance augmentation should be used, \code{"pred"} for posterior
#' predictive or \code{"real"} for realized (just the posterior).
#' @param ci.prob Probability for HPD credible intervals. Defaults to 0.95
#'
#' @return
#' A named list with names equal to the variables in \code{aggregation.data}. Each
#' element of the list is another list with elements:
#' \item{aggregated.abund}{The MCMC sample of the new aggregation}
#' \item{aggregation.summary}{A summary of the aggregation MCMC}
#' @export
#'
newAggregation <- function(fit, aggregation.data, type="pred", ci.prob=0.95){
  #require(coda)
  if(type == "pred") xxx <- fit$mcmc.sample$pred.site.abund
  if(type == "real") xxx <- fit$mcmc.sample$real.site.abund
  if(is.null(xxx)) stop("Site abundance data was not retained in the call to 'mcmc.aggreation()'\n Please re-run with 'keep.site.abund=TRUE'\n")
  site.name <- attr(fit,"site.name")
  time.name <- attr(fit,"time.name")
  site.idx <- data.frame(sapply(strsplit(colnames(xxx),"-"), function(x){paste(x[-1],collapse="-")}),
                         as.numeric(sapply(strsplit(colnames(xxx),"-"), function(x){x[[1]]})))
  colnames(site.idx) <- c(site.name, time.name)
  m1 <- merge(site.idx, aggregation.data, all=TRUE)
  ag.names <- colnames(aggregation.data)[colnames(aggregation.data)!=site.name]
  outlist <- vector("list", length(ag.names))
  #Tmat <- cbind(rep(1,length(unique(site.idx[,2]))), unique(site.idx[,2]))
  names(outlist) <- ag.names
  for(i in 1:length(ag.names)){
    cat("Processing", type, "aggregation:", ag.names[i], "...\n")
    a1 <- mcmc(t(apply(as.matrix(xxx), 1, FUN=function(v){aggregate(v, list(m1[,time.name],m1[,ag.names[i]]), FUN=sum)$x})))
    colnames(a1) <- apply(expand.grid(unique(site.idx[,time.name]), levels(factor(aggregation.data[,ag.names[i]]))),1,paste, collapse="-")
    ag.summary <- expand.grid(unique(site.idx[,time.name]), levels(factor(aggregation.data[,ag.names[i]])))
    colnames(ag.summary) <- c(time.name, ag.names[i])
    ag.summary$post.median.abund <- apply(a1,2,median)
    hpd <- HPDinterval(a1, ci.prob)
    ag.summary$low.hpd <- hpd[,1]
    ag.summary$hi.hpd <- hpd[,2]
    if(type=="pred") {
      outlist[[i]] <- list(mcmc.sample=list(aggregated.pred.abund=a1), aggregation.pred.summary=ag.summary)
    } else  {
      outlist[[i]] <- list(mcmc.sample=list(aggregated.real.abund=a1), aggregation.real.summary=ag.summary)
    }
    attr(outlist[[i]], "ln.adj") <- attr(fit, "ln.adj")
  }
  return(outlist)
}

