<!-- README.md is generated from README.Rmd. Please edit that file -->
[![DOI](https://zenodo.org/badge/8592075.svg)](https://zenodo.org/badge/latestdoi/8592075)

Fit regional trends to site-specific abundence data
---------------------------------------------------

This package fits a log-linear trend models to regions aggregated over
sites. The sites may contain missing surveys that are not temporally
aligned with the missing data at other sites, making direct aggregation
impossible. The functions within the package model the indivdual sites
with a semi-parametric (possibly, zero-inflated) model to interpolate
missing data from which regional aggregations can be made. By using
Markov Chain Monte Carlo, on can sample from the posterior predictive
distribution of the regional aggregations Then calculate the log-linear
trend over the time period of interest as a derived parameter. Using the
posterior predictive distribution allows incorporation of both parameter
uncertainty as well as uncertainty due to sampling the local abundance
processes.

### Disclaimer

*This software package is developed and maintained by scientists at the
NOAA Fisheries Alaska Fisheries Science Center and should be considered
a fundamental research communication. The reccomendations and
conclusions presented here are those of the authors and this software
should not be construed as official communication by NMFS, NOAA, or the
U.S. Dept. of Commerce. In addition, reference to trade names does not
imply endorsement by the National Marine Fisheries Service, NOAA. While
the best efforts have been made to insure the highest quality, tools
such as this are under constant development and are subject to change.*

### Example

Load packages for this example

``` r
library(agTrend.gam)
#> Loading required package: coda
#> Loading required package: mgcv
#> Loading required package: nlme
#> This is mgcv 1.8-25. For overview type 'help("mgcv-package")'.
#> Loading required package: tidyverse
#> ── Attaching packages ───────────────────────── tidyverse 1.2.1 ──
#> ✔ ggplot2 3.1.0     ✔ purrr   0.2.5
#> ✔ tibble  1.4.2     ✔ dplyr   0.7.7
#> ✔ tidyr   0.8.2     ✔ stringr 1.3.1
#> ✔ readr   1.1.1     ✔ forcats 0.3.0
#> ── Conflicts ──────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::collapse() masks nlme::collapse()
#> ✖ dplyr::filter()   masks stats::filter()
#> ✖ dplyr::lag()      masks stats::lag()
#> Loading required package: mvtnorm
#> Loading required package: furrr
#> Loading required package: future
#> 
#>  agTrend.gam 0.01.9000 (2018-10-24) 
#>  A demo is available at https://github.com/NMML/agTrend
```

Now we can load the data that is included with the `agTrend.gam` package
and filter it to the data we want for this example (i.e., 1990-2016).
Now we’ll add a photo method covariate to data (oblique photos prior to
2004 surveys = 1)

``` r
data(wdpsNonpups)
wdpsNonpups = wdpsNonpups %>% filter(YEAR>=1990) %>% mutate(obl = as.integer(YEAR<2004))
```

The data is then converted to the form necessary for the site abundance
imputing function

``` r
fit_data = wdpsNonpups %>% make.agtrend.data(abundance.name="COUNT", site.name="SITE", time.name="YEAR") 
```

Now, we count the number of positive counts at each site so that we can
remove sites that had only 1 positive count

``` r
fit_data <- fit_data %>% filter(num_nonzero>1)
```

The next step involves creating a prior distribution list for MCMC site
updating. An informative prior for the survey methodology correction is
obtained from analysis of another data set.

``` r
data("photoCorrection")
photoCorrection %>% mutate(log_ratio = log(OBLIQUE/VERTICAL)) -> photoCorrection
gamma_0 = photoCorrection %>% summarize(mean(log_ratio)) %>% as.double()
gamma_se = photoCorrection %>% summarize(sd(log_ratio)/sqrt(n())) %>% as.double()

obs.prior=list(gamma.mean=gamma_0, gamma.vcov=gamma_se^2)
```

For the linear trend at each site, a ridge regression penalty is used
such that the linear trend will be shrunk back to zero if the data do
not support a trend.

Now we begin the MCMC sampling using the `site.simulate` function. This
function performs the site augmentation and samples from the posterior
predictive distribution of the count data.

``` r
set.seed(123) 
fit_data = site.simulate(fit_data %>% slice(1:5), 
                         obs.formula=~obl-1,
                         min.k = 3, obs.prior=obs.prior, 
                         fit.only=F)
#> Fitting models to each site ...
#> Simulating missing observations ...
```

Now we add the regional designations to `fit_data` in order to summarize
the counts by desired regions

``` r
fit_data = fit_data %>% left_join(
  wdpsNonpups %>% select(SITE, REGION, RCA) %>% distinct(), by=c("site"="SITE")
  ) %>% mutate(TOTAL="TOTAL")
```

To demonstrate we choose to provide trends for the site region
designations. First, the counts are aggregated by region:

``` r
region_data = fit_data %>% ag.abundance(REGION)
```

The the abunance can be summarized with the following function call

``` r
region_summ = region_data %>% ag.summ()
```

Which can then be visualized with with the `ggplot2` package

``` r
p1 = ggplot(data = region_summ %>% filter(type=="prediction")) +
  geom_path(aes(y=estimate, x=time+1990)) + 
  geom_ribbon(
    aes(ymin=ci.lower, ymax=ci.upper,x=time+1990),alpha=0.2,fill="red"
    ) + facet_wrap(~REGION, ncol=2) +
  geom_pointrange(
    aes(x=time+1990, y=estimate, ymin=ci.lower, ymax=ci.upper),
    data=region_summ %>% filter(type=="realized")
    )
```

#### Disclaimer

<sub>This repository is a scientific product and is not official
communication of the Alaska Fisheries Science Center, the National
Oceanic and Atmospheric Administration, or the United States Department
of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project
code is provided on an ‘as is’ basis and the user assumes responsibility
for its use. AFSC-MML has relinquished control of the information and no
longer has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this
GitHub project will be governed by all applicable Federal law. Any
reference to specific commercial products, processes, or services by
service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation or favoring by the Department
of Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.</sub>
