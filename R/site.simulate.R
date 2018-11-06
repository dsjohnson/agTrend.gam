#' @title Posterior predictive sampling, aggregtion of abundance counts, and linear trend summary
#'
#' @description Function for sampling from the posterior predictive distribution of abundance (counts) at
#' individual sites. Then aggregating the counts over the specified aggregation variable.
#'
#' @param data
#' @param obs.formula
#' @param sig.abund
#' @param min.k
#' @param forecast
#' @param family
#' @param ln.adj
#' @param prior.list
#' @param upper
#' @param lower
#' @param fit.only
#' @param post.size
#' @param run.parallel
#'
#' @return
#' A named list with the following elements:
#' \item{site.summary}{A summary of the abundance augmentation for every site and every time.}
#' \item{posterior.sample}{A named list containing all of the MCMC sample after thinning.}
#' \item{original.data}{The original data in \code{data}.}
#'
#' @export
#' @import truncnorm
#' @import mvtnorm
#' @import mgcv
#'
site.simulate <- function(data, obs.formula=NULL, sig.abund=NULL,
                          min.k = 3, forecast=0,
                          family=c("tweedie","log.normal"),
                          ln.adj=0,
                          obs.prior=NULL,
                          upper=Inf, lower=-Inf,
                          fit.only=FALSE,
                          post.size=1000,
                          run.parallel=FALSE
){


  if(run.parallel){
    plan(multiprocess)
  } else plan(sequential)

  message("Fitting models to each site ... ")
  data <- data %>%
    mutate(
      attempt = furrr::future_map(data, site.fit.gam %>% safely(),
                                  obs.formula=obs.formula, min.k=min.k,
                                  ln.adj=ln.adj, obs.prior=obs.prior,
                                  .progress = TRUE)
    )
  data <- data %>% mutate(
    fits = purrr::map(attempt, ~.x[["result"]]),
    fit.error = purrr::map_chr(attempt, ~ifelse(!is.null(.x[["error"]]$message), .x[["error"]]$message, NA_character_))
  ) %>% select(-attempt)
  if(fit.only){
    return(data)
  } else {
    if(!all(is.na(data$fit.error))){
      stop("There were model fitting errors! run again with 'fit.only=TRUE' to diagnose.")
    } else{
      message("Simulating missing observations ... ")
      max.time = purrr::map_dbl(data$data, ~{max(.x$time)}) %>% max()
      data = data %>% mutate(
        post = furrr::future_map2(fits, data,
                                  site.sim.gam %>% safely(), max.time=max.time+forecast,
                                  size=post.size,
                                  .progress = TRUE)
      )
      data <- data %>% mutate(
        post.sample = purrr::map(post, ~.x[["result"]]),
        post.error = purrr::map_chr(post, ~ifelse(!is.null(.x[["error"]]$message), .x[["error"]]$message, NA_character_))
      ) %>% select(-post)
      data<-data %>%
        mutate(
          post.sample.real = purrr::map2(data, post.sample, sim.to.realized %>% possibly("NA"))
        )
      return(data)
    }
  }
  plan(sequential)
}

