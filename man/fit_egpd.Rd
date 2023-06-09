% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idf_functions.R
\name{fit_egpd}
\alias{fit_egpd}
\title{Function to fit the three parameter EGPD of Naveau et al 2016}
\usage{
fit_egpd(
  data,
  fitting_method = "mle",
  init,
  declustering_duration = 1,
  init_time_step = 1,
  left_censoring_value = 0,
  plot = FALSE,
  optim_algo = "Nelder-Mead"
)
}
\arguments{
\item{data}{vector of precipitation intensities}

\item{fitting_method}{eiither \code{"mle"} for maximum likelihood or \code{"pwm"} for probability weighted momoments. Defaults to \code{"mle"}}

\item{init}{a vector of length three, respectively the initial values for \code{kappa, sigma, xi}. Can be omitted.}

\item{declustering_duration}{a scalar .Whether the data should be temporally declustered, if yes the time step for each duration, default to a value of  1 otherwise}

\item{init_time_step}{a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd timestep will be selected, and then a sequence is applied}

\item{left_censoring_value}{a scalar or vector of length(durations), defaults to \code{0}. The left censoring value to be applied to data of each duration. If a scalar is given, it will be divided (scaled) by durations. It will be ignored if  \code{auto_fit = TRUE}.}

\item{plot}{logical, defaults to False, if \code{TRUE}, qqplot will be displayed, and a ggplot object will be returned.}

\item{optim_algo}{the \code{optim} algorithm to use. defaults to \code{"Nelder-Mead" }}
}
\value{
# A list:
   fitted parameters: the three egpd paramters
   likelihood: the likelihood values
   AIC: the AIC of the model
   nRMSE:  root-mean-square-error, normalized by the mean of \code{data}, to assess the quality of the fit. The smaller the better, a value of <= 0.1 is always good
   qqplot: a ggplot object, the quantile quantile plot.
   inference: a charcter, either "mle" or "pwm", for the fitting method used
}
\description{
The inference is done with either maximum likelihood (mle) or probability weighted method (pwm).
}
\details{
to be added
}
\examples{
 ## load the data
 data("precipdata")

\dontrun{
## fit EGPD to the data, check the qqplot and the nRMSE
fit_egpd(precipdata$SCH, declustering_duration = 3,  left_censoring_value = 0 , plot = T)

## try left censoring
fit_egpd(precipdata$SCH, declustering_duration = 3,  left_censoring_value = 1 , plot = T)

## try a couple of left censoring values, and check the nRMSE, and the corresponding qqplot
lapply(seq(1,3,0.5), function(x){
  fit = fit_egpd(precipdata$SCH, declustering_duration = 3,  left_censoring_value = x ,
  plot = T)
  print(fit$nRMSE)
  fit$qqplot +
   labs(subtitle = paste0("nRMSE = ", round(fit$nRMSE, 3)))
})

## Now try fiitng using PWM method
fit_egpd(precipdata$SCH, fitting_method = "pwm", declustering_duration = 3,
left_censoring_value = 0 , plot = T)
}

}
