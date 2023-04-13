#' Intensity Duration Frequency model (IDF) using the Extended Generalised Pareto Distribution
#'
#'
#' IDF curves are usually modelled using Generalized Exteme Value distributions (GEV).
#' Here, we consider the three parameter EGPD of Naveau et al 2016 to bulid IDF curves. The package provides the implemenattion
#' of 10 diffent IDF modelling approaches as contained in Haruna et al 2023.
#'
#' @docType package
#'
#' @author Abubakar HARUNA \email{abbaharuna87@gmail.com}
#'
#' @name egpdIDF
#'
#' @import mev
#' @import segmented
#' @import doParallel
#' @import parallel
#' @importFrom extRemes fevd
#' @import ggplot2
#' @importFrom scales trans_breaks trans_format math_format
#' @import  stats
#' @importFrom dplyr group_by mutate summarise_at select rename  summarise arrange group_split
#' @import  purrr
#' @import  RcppRoll
#' @import  tidyr
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom graphics legend lines
NULL

globalVariables(c("duration", "i", "x", "y", 'q.mle.U', 'q.mle.L','tr_years', 'tr_quantiles', '.x', 'dur', 'tr_idf', 'Tr'))

#'A fictious hourly  time series of precipitation
#'
#' @format ## `precipdata`
#' A data frame of two columns and 1577952 rows:
#' \describe{
#'   \item{date}{date}
#'   \item{SCH}{precipitation amounts}
#' }
"precipdata"
