#' @title Calculate log-linear trend estmates from aggregated abundance data
#' @description Takes aggregated abundance values calculated from `ag.abundance` and
#' calculates trend estimates for each region
#' @param object A tibble produced by `ag.abundance()`
#' @param time.range A length 2 vector giving the start and end times for the
#' trend estimates. These times are between 0 and the number times in the data set.
#' Not, for example, specific years.
#' @import dplyr modeest
#' @author Devin S. Johnson
#' @export

ag.trend = function(object, time.range=NULL, ci.prob=0.95){
  get_mode = function(x){modeest::mlv(x, method="shorth")$M}
  if(is.null(time.range)) time.range=c(0,nrow(object$post.sample[[1]])-1)
  times = time.range[1]:time.range[2]
  trends = object %>% select(post.sample) %>%
    mutate(
      post.sample = map(post.sample, ~{.x[times+1,]}),
      post.list = map(post.sample, ~{split(.x, rep(1:ncol(.x), each = nrow(.x)))})
    ) %>% select(-post.sample)
  trends = trends %>% mutate(
    post.list = purrr::map(post.list,
                           ~{
                             purrr::map(.x, ~{data.frame(y=.x, trend=times)})
                           })
  )
  trends = trends %>%
    mutate(
      trend.sample = purrr::map(post.list,
                                ~{
                                  map(.x, ~{glm(y~trend, family="poisson", data=.x)$coefficients}) %>%
                                    simplify2array() %>% t()
                                }
      )
    ) %>% select(-post.list) %>% bind_cols(object[,1],.)
  trends = trends  %>%
    mutate(
      growth.estimate = map_dbl(trend.sample, ~{get_mode(100*(exp(.x[,2])-1))}),
      growth.ci = map(trend.sample,
                      ~{
                        100*(exp(.x[,2])-1) %>% coda::mcmc() %>%
                          coda::HPDinterval() %>% as.data.frame
                      }
      )) %>% unnest(growth.ci) %>%
    rename(growth.lower=lower, growth.upper=upper)
  trends = trends %>% mutate(
    trend.line = map(trend.sample, ~{
      tl.sample=cbind(1,times)%*%t(.x) %>% exp()
      fitted = apply(tl.sample,1,get_mode)
      out=data.frame(time=times, fitted)
      trend.line.ci = tl.sample %>% t() %>%
        coda::mcmc() %>% coda::HPDinterval() %>% as.data.frame()
      colnames(trend.line.ci)=paste0("trend.line.",colnames(trend.line.ci))
      cbind(out, trend.line.ci)
    }
    )
  )
  return(trends)
}
