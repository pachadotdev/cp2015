#' @title Caliendo and Parro (2015) Quantitative Trade Model
#' @useDynLib cp2015, .registration = TRUE
#' @importFrom data.table ':=' as.data.table copy data.table set setnames setorderv setDT
#' @keywords internal
"_PACKAGE"

#' @title Data to Replicate Section 5.1
#' @description To quantify the trade and welfare effects of NAFTA. Welfare
#' effects from NAFTA's tariff reductions while fixing the tariff to and from
#' the rest of the world to the year 1993.
#'
#' @format ## `cp2015nafta`
#' A list of seven elements. Sets: regions (e.g., Brazil) and sectors (e.g.,
#' Agriculture). Intermediate consumption: input (e.g., Agriculture), sector
#' (e.g, Food), region, and value. Final consumption: sector, region, and value.
#' Value added: Sector, region, and value. Trade: sector, exporter, importer,
#' value, tariff, tariff baseline, tariff counterfactual, and relative changes
#' of the iceberg trade costs. Deficit: region and value. Theta (e.g., trade
#' elasticities): sector and value.
#'
#' @source Estimates of the Trade and Welfare Effects of NAFTA (Caliendo and Parro, 2015)
"cp2015nafta"
