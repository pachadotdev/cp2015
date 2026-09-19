
# cp2015

The goal of `cp2015` package is to implement the Caliendo and Parro (2015) quantitative trade model in R.

## Installation

To install the development version, run the following command:

``` r
remotes::install_github("pachadotdev/cp2015")
```

### Backend

This port uses `data.table` for data preparation and result aggregation, and `armadillo4r`/`cpp4r` for the numerical
solver. It no longer depends on `dplyr`, `magrittr`, `Rcpp`, or `xtensor`.

## Example

This replicates the simulation presented in the subsection 5.1 of the article
["Estimates of the Trade and Welfare Effects of NAFTA"](https://academic.oup.com/restud/article/82/1/1/1547758).

```r
library(cp2015)
library(data.table)

options(knitr.kable.NA = "")
```

## Data

The package includes `cp2015nafta`, a list containing the required model sets,
intermediate consumption, final consumption, value added, bilateral trade and
tariff data, deficits, and trade elasticities.

```r
str(cp2015nafta)
```

```r
List of 7
 $ sets                    :List of 2
  ..$ regions: chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ sectors: chr [1:40] "Agriculture" "Mining" "Food" "Textile" ...
 $ intermediate_consumption:Classes ‘data.table’ and 'data.frame':      49600 obs. of  4 variables:
  ..$ input : chr [1:49600] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ sector: chr [1:49600] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ region: chr [1:49600] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ value : num [1:49600] 2.62e+09 2.96e+09 5.20e+08 3.89e+09 7.63e+09 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)> 
 $ final_consumption       :Classes ‘data.table’ and 'data.frame':      1240 obs. of  3 variables:
  ..$ sector: chr [1:1240] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ region: chr [1:1240] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ value : num [1:1240] 7.16e+09 6.48e+09 4.76e+09 1.59e+10 4.74e+09 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)> 
 $ value_added             :Classes ‘data.table’ and 'data.frame':      1240 obs. of  3 variables:
  ..$ sector: chr [1:1240] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ region: chr [1:1240] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ value : num [1:1240] 1.24e+10 1.18e+10 5.12e+09 2.83e+10 1.51e+10 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)> 
 $ trade                   :Classes ‘data.table’ and 'data.frame':      38440 obs. of  9 variables:
  ..$ sector    : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ exporter  : chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ importer  : chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ value     : num [1:38440] 1.91e+10 2.19e+06 5.56e+04 8.87e+07 1.50e+06 ...
  ..$ tariff    : num [1:38440] 0 0.0417 0.0408 0.0488 0.0288 ...
  ..$ tariff_bln: num [1:38440] 0 0.0417 0.0408 0.0488 0.0288 ...
  ..$ tariff_cfl: num [1:38440] 0 0.0417 0.0408 0.0488 0.0288 ...
  ..$ d_bln     : num [1:38440] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ d_cfl     : num [1:38440] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)> 
 $ deficit                 :Classes ‘data.table’ and 'data.frame':      31 obs. of  2 variables:
  ..$ region: chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ D     : num [1:31] 3.51e+09 -4.19e+09 1.62e+10 -1.43e+10 -1.07e+10 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)> 
 $ theta                   :Classes ‘data.table’ and 'data.frame':      40 obs. of  2 variables:
  ..$ sector: chr [1:40] "Agriculture" "Mining" "Food" "Textile" ...
  ..$ value : num [1:40] 9.11 13.53 2.62 8.1 11.5 ...
  ..- attr(*, ".internal.selfref")=<pointer: (nil)
```

## Welfare effects from NAFTA's tariff reductions (table 2)

Caliendo and Parro (2015) impose zero aggregate deficits in the baseline and
counterfactual scenarios. Use `zero_aggregate_deficit = TRUE` to reproduce that
scenario.

```r
results <- run_cp2015(data = cp2015nafta, zero_aggregate_deficit = TRUE)
str(results)
```

```r
Solving the baseline.
Iteration: 1 - ||F(x)||: 4075219334
Iteration: 100 - ||F(x)||: 0.009735903192
Iteration: 200 - ||F(x)||: 0.001738052373
Iteration: 300 - ||F(x)||: 0.0003572560362
Iteration: 400 - ||F(x)||: 6.213691551e-05
Iteration: 500 - ||F(x)||: 1.315719258e-05
Iteration: 600 - ||F(x)||: 2.223333097e-06
Iteration: 700 - ||F(x)||: 4.8383701e-07
Iteration: 796 - ||F(x)||: 9.686823042e-08
Solving the counterfactual.
Iteration: 1 - ||F(x)||: 4074074133
Iteration: 100 - ||F(x)||: 0.01007323017
Iteration: 200 - ||F(x)||: 0.001681247597
Iteration: 300 - ||F(x)||: 0.0003469702056
Iteration: 400 - ||F(x)||: 5.838458239e-05
Iteration: 500 - ||F(x)||: 1.212108083e-05
Iteration: 600 - ||F(x)||: 2.028661454e-06
Iteration: 700 - ||F(x)||: 4.233016861e-07
Iteration: 779 - ||F(x)||: 9.651468184e-08

List of 13
 $ c_nj_hat        :Classes ‘data.table’ and 'data.frame':      1240 obs. of  3 variables:
  ..$ region: chr [1:1240] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ sector: chr [1:1240] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ c_hat : num [1:1240] 0.999 0.999 0.999 0.999 0.998 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ P_nj_hat        :Classes ‘data.table’ and 'data.frame':      1240 obs. of  3 variables:
  ..$ region: chr [1:1240] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ sector: chr [1:1240] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ P_hat : num [1:1240] 0.999 0.999 0.999 0.999 0.998 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ pi_nij          :Classes ‘data.table’ and 'data.frame':      38440 obs. of  5 variables:
  ..$ importer: chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ exporter: chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ sector  : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ pi_bln  : num [1:38440] 9.85e-01 6.35e-05 5.84e-04 1.43e-02 8.96e-04 ...
  ..$ pi_cfl  : num [1:38440] 9.85e-01 6.37e-05 5.82e-04 1.43e-02 8.88e-04 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ X_nj            :Classes ‘data.table’ and 'data.frame':      1240 obs. of  4 variables:
  ..$ region: chr [1:1240] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ sector: chr [1:1240] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ X_bln : num [1:1240] 1.89e+10 1.78e+10 1.06e+10 4.10e+10 3.28e+10 ...
  ..$ X_cfl : num [1:1240] 1.88e+10 1.78e+10 1.06e+10 4.10e+10 3.32e+10 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ I_n             :Classes ‘data.table’ and 'data.frame':      31 obs. of  3 variables:
  ..$ region: chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ I_bln : num [1:31] 2.10e+11 3.46e+11 1.75e+11 4.62e+11 5.72e+11 ...
  ..$ I_cfl : num [1:31] 2.10e+11 3.46e+11 1.74e+11 4.61e+11 5.69e+11 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ P_n_hat         :Classes ‘data.table’ and 'data.frame':      31 obs. of  2 variables:
  ..$ region : chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ P_n_hat: num [1:31] 0.999 0.999 0.999 0.999 0.996 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ w_n_hat         :Classes ‘data.table’ and 'data.frame':      31 obs. of  2 variables:
  ..$ region: chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ w_hat : num [1:31] 0.999 0.999 0.999 0.999 0.999 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ trade           :Classes ‘data.table’ and 'data.frame':      38440 obs. of  9 variables:
  ..$ importer : chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ exporter : chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ sector   : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ tau_bln  : num [1:38440] 0 0.0371 0.0522 0.0941 0.000467 ...
  ..$ tau_cfl  : num [1:38440] 0 0.0371 0.0522 0.0941 0.000467 ...
  ..$ d_bln    : num [1:38440] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ d_cfl    : num [1:38440] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ trade_bln: num [1:38440] 1.86e+10 1.09e+06 5.90e+06 5.37e+08 2.94e+07 ...
  ..$ trade_cfl: num [1:38440] 1.86e+10 1.10e+06 5.87e+06 5.36e+08 2.94e+07 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ tot             :Classes ‘data.table’ and 'data.frame':      38440 obs. of  4 variables:
  ..$ partner: chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ region : chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ sector : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ tot    : num [1:38440] 0.00 5.86e-08 -2.83e-06 -2.31e-04 -1.29e-05 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ vot             :Classes ‘data.table’ and 'data.frame':      38440 obs. of  4 variables:
  ..$ region : chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ partner: chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ sector : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ vot    : num [1:38440] 0.00 3.13e-08 -5.72e-07 -2.35e-06 3.83e-09 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ tech            :Classes ‘data.table’ and 'data.frame':      38440 obs. of  4 variables:
  ..$ region : chr [1:38440] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ partner: chr [1:38440] "Argentina" "Argentina" "Argentina" "Argentina" ...
  ..$ sector : chr [1:38440] "Agriculture" "Agriculture" "Agriculture" "Agriculture" ...
  ..$ tech   : num [1:38440] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
 $ welfare         :Classes ‘data.table’ and 'data.frame':      31 obs. of  6 variables:
  ..$ region  : chr [1:31] "Argentina" "Australia" "Austria" "Brazil" ...
  ..$ total   : num [1:31] 0.000978 0.000159 -0.004278 -0.002096 -0.063816 ...
  ..$ tot     : num [1:31] 0.000254 0.000385 -0.002212 -0.00218 -0.108102 ...
  ..$ vot     : num [1:31] 7.24e-04 -2.26e-04 -2.07e-03 8.38e-05 4.43e-02 ...
  ..$ realwage: num [1:31] 0.001313 0.000662 -0.002147 -0.001851 0.32283 ...
  ..$ tech    : num [1:31] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, ".internal.selfref")=<pointer: 0x5559f29a1f00>
 $ convergence_info:Classes ‘data.table’ and 'data.frame':      2 obs. of  3 variables:
  ..$ scenario      : chr [1:2] "Baseline" "Counterfactual"
  ..$ criteria_value: num [1:2] 9.69e-08 9.65e-08
  ..$ message       : chr [1:2] "Successful convergence" "Successful convergence"
  ..- attr(*, ".internal.selfref")=<pointer: 0x5565cc0d4f00> 
```

The welfare decomposition includes terms of trade, volume of trade, and
technical efficiency. The latter is zero in this tariff-only example.

```r
nafta <- c("Mexico", "Canada", "USA")
welfare <- results$welfare[region %in% nafta]
welfare[, tech := NULL]
welfare <- welfare[match(nafta, region), ]

knitr::kable(
  welfare,
  digits = 2,
  col.names = c("Country", "Total", "Terms of trade", "Volume of trade", "Real wages"),
  caption = "Welfare effects from NAFTA's tariff reductions"
)
```

```
Table: Welfare effects from NAFTA's tariff reductions

|Country | Total| Terms of trade| Volume of trade| Real wages|
|:-------|-----:|--------------:|---------------:|----------:|
|Mexico  |  1.31|          -0.41|            1.72|       1.72|
|Canada  | -0.06|          -0.11|            0.04|       0.32|
|USA     |  0.08|           0.04|            0.04|       0.11|
```

## Bilateral welfare effects from NAFTA's tariff reductions (table 3)

```r
tot_bilateral <- results$tot[region %in% nafta]
tot_bilateral[, partner := fifelse(partner %in% nafta, "NAFTA", "ROW")]
tot_bilateral <- tot_bilateral[, .(ToT = sum(tot)), by = .(region, partner)]

vot_bilateral <- results$vot[region %in% nafta]
vot_bilateral[, partner := fifelse(partner %in% nafta, "NAFTA", "ROW")]
vot_bilateral <- vot_bilateral[, .(VoT = sum(vot)), by = .(region, partner)]

bilateral <- merge(
  tot_bilateral,
  vot_bilateral,
  by = c("region", "partner"),
  all = TRUE,
  sort = FALSE
)

bilateral <- dcast(
  bilateral,
  region ~ partner,
  value.var = c("ToT", "VoT")
)

bilateral <- bilateral[match(nafta, region), ]

knitr::kable(
  bilateral,
  digits = 2,
  col.names = c("Country", "ToT NAFTA", "ToT ROW", "VoT NAFTA", "VoT ROW"),
  caption = "Bilateral welfare effects from NAFTA's tariff reductions"
)
```

```
Table: Bilateral welfare effects from NAFTA's tariff reductions

|Country | ToT NAFTA| ToT ROW| VoT NAFTA| VoT ROW|
|:-------|---------:|-------:|---------:|-------:|
|Mexico  |     -0.39|   -0.02|      1.80|   -0.08|
|Canada  |     -0.09|   -0.02|      0.08|   -0.04|
|USA     |      0.03|    0.01|      0.04|    0.00|
```

## Sectoral contribution to welfare effects from NAFTA's tariff reductions (table 4)

```r
tot_sector <- results$tot[region %in% nafta]
tot_sector <- tot_sector[, .(ToT = sum(tot)), by = .(region, sector)]
tot_sector[, ToT := 100 * ToT / sum(ToT), by = region]

vot_sector <- results$vot[region %in% nafta]
vot_sector <- vot_sector[, .(VoT = sum(vot)), by = .(region, sector)]
vot_sector[, VoT := 100 * VoT / sum(VoT), by = region]

sectoral <- merge(
  tot_sector,
  vot_sector,
  by = c("region", "sector"),
  all = TRUE,
  sort = FALSE
)

sectoral <- dcast(
  sectoral,
  sector ~ region,
  value.var = c("ToT", "VoT")
)

setcolorder(
  sectoral,
  c(
    "sector",
    "ToT_Mexico", "VoT_Mexico",
    "ToT_Canada", "VoT_Canada",
    "ToT_USA", "VoT_USA"
  )
)

knitr::kable(
  sectoral,
  digits = 2,
  col.names = c("Sector", gsub("_", " ", colnames(sectoral)[-1])),
  caption = "Sectoral contribution to welfare effects from NAFTA's tariff reductions"
)
```

```
Table: Sectoral contribution to welfare effects from NAFTA's tariff reductions

|Sector          | ToT Mexico| VoT Mexico| ToT Canada| VoT Canada| ToT USA| VoT USA|
|:---------------|----------:|----------:|----------:|----------:|-------:|-------:|
|Agriculture     |      -0.13|       2.87|       3.41|      -0.01|    3.41|    0.65|
|Air Transport   |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Auto            |      13.79|       4.78|      29.50|      27.77|   15.85|    4.47|
|Aux Transport   |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Basic metals    |       1.07|       3.02|      10.08|       1.48|    3.40|    1.05|
|Chemicals       |       0.57|       2.15|       5.74|       0.08|    5.60|    1.11|
|Communication   |      20.97|       3.64|       2.67|       0.15|   11.61|    4.58|
|Computer        |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Construction    |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Education       |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Electrical      |      41.20|      25.77|       1.37|       7.18|   24.19|   42.22|
|Electricity     |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Finance         |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Food            |       0.45|       1.17|       3.56|       2.37|    3.16|    1.04|
|Health          |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Hotels          |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Land Transport  |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Machinery n.e.c |       3.68|       4.32|       5.16|      -0.02|    5.63|    0.65|
|Medical         |       4.72|       1.34|       0.94|      -0.23|    3.48|    4.46|
|Metal products  |       0.90|       5.56|       2.22|       7.99|    1.61|    1.06|
|Minerals        |       0.05|       0.73|       0.93|       0.47|    0.70|    0.57|
|Mining          |      -3.01|       0.25|       4.04|      -0.20|    1.54|    0.04|
|Office          |       8.37|       4.72|       2.32|      -0.82|    3.50|    1.43|
|Other           |       2.63|       1.92|       0.81|      -0.10|    2.90|    1.69|
|Other Business  |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Other Transport |       0.21|       0.82|      12.93|      -0.97|    1.51|    0.32|
|Other services  |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Paper           |       0.39|       3.82|       5.86|       0.49|    2.83|    0.33|
|Petroleum       |      -0.09|      14.61|       0.60|      30.41|    1.85|   11.37|
|Plastic         |       0.62|       4.21|       2.53|       7.56|    1.61|    0.32|
|Post            |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Private         |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Public          |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|R&D             |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Real State      |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Renting Mach    |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Retail          |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Textile         |       3.30|      12.03|       1.15|      16.19|    4.32|   22.24|
|Water Transport |       0.00|       0.00|       0.00|       0.00|    0.00|    0.00|
|Wood            |       0.30|       2.26|       4.17|       0.24|    1.31|    0.41|
```

## Trade effects from NAFTA's tariff reductions (table 5)

```r
trade <- results$trade
trade <- trade[importer %in% nafta & exporter %in% nafta]
trade <- trade[, .(trade_bln = sum(trade_bln), trade_cfl = sum(trade_cfl)),
  by = .(importer, exporter)
]
trade[, change := (trade_cfl / trade_bln - 1) * 100]
trade[exporter == importer, change := NA_real_]

trade <- dcast(
  trade,
  importer ~ exporter,
  value.var = "change"
)

trade <- trade[match(nafta, importer), ]
setcolorder(trade, c("importer", "Mexico", "Canada", "USA"))

knitr::kable(
  trade,
  digits = 2,
  col.names = c("Importer/Exporter", nafta),
  caption = "Trade effects from NAFTA's tariff reductions"
)
```

```
Table: Trade effects from NAFTA's tariff reductions

|Importer/Exporter | Mexico| Canada|    USA|
|:-----------------|------:|------:|------:|
|Mexico            |       | 116.60| 118.31|
|Canada            |  58.57|       |   9.49|
|USA               | 109.54|   6.57|       |
```
