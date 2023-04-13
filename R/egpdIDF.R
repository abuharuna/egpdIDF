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
#' @import  dplyr
#' @import  stats
#' @import  purrr
#' @import  RcppRoll
#' @import  tidyr
#' @import foreach
#' @importFrom graphics legend lines
NULL

globalVariables(c("duration", "i"))
