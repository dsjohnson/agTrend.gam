#' @title Produce data object for use in \code{site.simulate}
#'
#' @description This function creates the specific data form used my the remaining functions in the
#' package for simulating missing information at all sites.
#'
#' @param object
#' @param abundance.name
#' @param site.name
#' @param time.name
#'
#' @return A tibble
#' @author Devin S. Johnson
#'
#' @import dplyr rlang
#' @export

make.agtrend.data <- function(object, abundance.name, site.name, time.name){
  abundance.name = rlang::enquo(abundance.name)
  site.name = rlang::enquo(site.name)
  time.name = rlang::enquo(time.name)
  object = object %>% filter(!is.na(!!abundance.name))
  object %>% dplyr::mutate(
    y := !! abundance.name,
    site := !! site.name,
    time := !! time.name
    ) %>%
    dplyr::mutate(time=time-min(time),
                  ym1 = ifelse(y>0, y-1, NA),
                  z = 1.0*(y>0)
    ) %>% dplyr::group_by(site) %>% tidyr::nest() -> object
  object %>% mutate(
    n_survey = purrr::map_int(data, nrow),
    num_nonzero = purrr::map_dbl(data, ~{sum(.x$z)}),
    max_abund = purrr::map_dbl(data, ~{max(.x$y)}),
    avg_abund = purrr::map_dbl(data, ~{mean(.x$y, na.rm=TRUE)})
  ) -> object
  return(object)
}
