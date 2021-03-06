
#' R package for fitting temporal trends to abundence data aggregated over large regions when subregions have missing data
#'
#' This package fits a log-linear trend models to regions aggregated over sites. The sites may contain missing surveys that
#' are not temporally aligned with the missing data at other sites, making direct aggregation impossible. The functions within the package
#' model the indivdual sites with a semi-parametric (possibly, zero-inflated) model to interpolate missing data from which regional aggregations
#' can be made. By using Monte Carlo approach, on can sample from the posterior predictive distribution of the regional aggregations
#' Then calculate the log-linear trend over the time period of interest as a derived parameter. Using the posterior predictive distribution
#' allows incorporation of both parameter uncertainty as well as uncertainty due to sampling the local abundance processes.
#'
#' \tabular{ll}{ Package: \tab agTrend.gam\cr
#' Type: \tab Package\cr
#' Version: \tab 0.01.9000\cr
#' Date: \tab 2018-10-24\cr
#' License: \tab CC0\cr
#' LazyLoad: \tab
#' yes\cr }
#'
#' @name agTrend.gam-package
#' @aliases agTrend.gam-package agTrend.gam
#' @docType package
#' @author Devin S. Johnson
#'
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#'
#' @importFrom stats aggregate as.formula binomial coef dnorm glm lm median model.matrix
#' pnorm predict quantile residuals rgamma rnorm var
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

#' Steller sea lion survey data (nonpup counts) collected by the NOAA National Marine Mammal Laboratory in the western Distinct
#' Population Segment of the stock (wDPS).
#'
#'
#' @name wdpsNonpups
#' @docType data
#' @format
#'
#' A data frame with 4118 observations on the following 5 variables.
#'
#' \describe{
#' \item{SITE}{Site where survey was taken}
#'
#' \item{REGION}{Region in which the site is located}
#'
#' \item{RCA}{Rookery Cluster Area in which the site belongs}
#'
#' \item{YEAR}{Year in which the survey was conducted}
#'
#' \item{COUNT}{The count of nonpups observed}
#' }
#'
#' @references Need a tech report / memo reference...
#'
#' @source Alaska Ecosystems Program, Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(wdpsNonpups)
#' head(wdpsNonpups)
NULL

#' Data from a pilot study of collecting Steller sea lion counts from photos taken by hand at
#' oblique angles vs. vertical medium-format photos.
#'
#'
#' @name photoCorrection
#' @docType data
#' @format
#'
#' A data frame with 20 observations on the following 4 variables.
#'
#' \describe{
#' \item{SITE}{Site where survey was taken}
#'
#' \item{REGION}{Region in which the site is located}
#'
#' \item{OBLIQUE}{Survey count taken from oblique photo source}
#'
#' \item{VERTICAL}{Survey count taken from vertical medium-format photo source}
#'
#' }
#'
#' @source Alaska Ecosystems Program, Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(photoCorrection)
#' head(photoCorrection)
NULL

#' Pup counts from surveys of the eastern Distict Populations Segment (eDPS) of Steller sea lions.
#'
#' @name edpsNonpups
#' @docType data
#' @format
#'
#' A data frame with 86 observations on the following 3 variables.
#'
#' \describe{
#' \item{SITE}{Site where the survey was taken.}
#'
#' \item{REGION}{One of 4 eDPS subregions: southeast Alaska (SE AK), British Columbia (BC), Oregon (OR), and California (CA)}
#'
#' \item{YEAR}{Year the survey was conducted}
#'
#' \item{COUNT}{Aggregated count of pups of sites within each subregion.}
#' }
#'
#' @source Alaska Ecosystems Program, Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(edpsNonpups)
#' head(edpsNonpups)
NULL


#' Pup counts from surveys of the eastern Distict Populations Segment (eDPS) of Steller sea lions.
#'
#' @name edpsPups
#' @docType data
#' @format
#'
#' A data frame with 34 observations on the following 3 variables.
#'
#' \describe{
#' \item{SITE}{Site where the survey was conducted.}
#'
#' \item{REGION}{One of 4 eDPS subregions: southeast Alaska (SE AK), British Columbia (BC), Oregon (OR), and California (CA)}
#'
#' \item{YEAR}{Year the survey was conducted}
#'
#' \item{COUNT}{Aggregated count of pups of sites within each subregion.}
#' }
#'
#' @source Alaska Ecosystems Program, Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(edpsPups)
#' head(edpsPups)
NULL

#' Pup counts from surveys of the western Distict Populations Segment (wDPS) of Steller sea lions.
#'
#' @name wdpsPups
#' @docType data
#' @format
#'
#' A data frame with 694 observations on the following 5 variables.
#'
#' \describe{
#' \item{SITE}{Surveyed sites in the wDPS}
#'
#' \item{REGION}{Region of the surveyed site}
#'
#' \item{RCA}{Rookery cluster area of the surveyed site}
#'
#' \item{YEAR}{Survey year}
#'
#' \item{COUNT}{Survey count of pups at each site}
#' }
#'
#' @source Alaska Ecosystems Program, National Marine Mammal Laboratory, Alaska
#' Fisheries Science Center, National Marine Fisheries Service, NOAA, 7600 Sand
#' Point Way, NE Seattle, WA 98115
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(wdpsPups)
#' head(wdpsPups)
NULL



.onAttach <- function(library, pkgname)
  {
    info <-utils::packageDescription(pkgname)
    package <- info$Package
    version <- info$Version
    date <- info$Date
    packageStartupMessage(
      paste("\n",paste(package, version, paste("(",date, ")", sep=""), "\n"),
            "A demo is available at https://github.com/dsjohnson/agTrend.gam"
            )
    )
  }

