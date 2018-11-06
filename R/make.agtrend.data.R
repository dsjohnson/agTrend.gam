#' @title Produce data object for use in \code{site.simulate}
#'
#' @description This function creates the specific data form used my the remaining functions in the
#' package for simulating missing information at all sites.
#'
#' @param object
#' @param abundance.name
#' @param site.name
#' @param time.name
#' @param gam.cut
#'
#' @return A tibble
#' @author Devin S. Johnson
#'
#' @import dplyr
#' @export

make.agtrend.data <- function(object, abundance.name, site.name, time.name, gam.cut){
  object %>% dplyr::mutate(
    y=.data[[abundance.name]],
    site=.data[[site.name]],
    time=.data[[time.name]]) %>%
    dplyr::mutate(time=time-min(time),
                  ym1 = ifelse(y>0, y-1, NA),
                  z = 1.0*(y>0)
    ) %>% dplyr::group_by(site) %>% tidyr::nest() -> object
  object %>% mutate(
    n_survey = purrr::map_int(data, nrow),
    num_nonzero = purrr::map_dbl(data, ~{sum(.x$z)}),
    max_abund = purrr::map_dbl(data, ~{max(.x$y)})
  ) -> object
  return(object)
}
