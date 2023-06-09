% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{aggregate_data}
\alias{aggregate_data}
\title{function to aggregate data into different intensities of longer durations}
\usage{
aggregate_data(sample_data, st_code, durations)
}
\arguments{
\item{sample_data}{data_frame of at least two columns, one of which is  named \code{date} containg the dates, should be in correct \code{R} date format. The other column contains the timeseries for each station}

\item{st_code}{a character. Must be one of \code{colnames}  of \code{sample_data}.}

\item{durations}{a numerical vector of aggregation durations. Must be in hours, if the resolution is hourly or in minutes. If in days, it should be in days, eg 1,2,3}
}
\value{
# A data frame: with \code{ncol = length(durations)}. Each containing the aggregated time series
}
\description{
Given a precipitation data of particular duration, the function returns a dataframe containing the aggregate data
}
\details{
to be added
}
\examples{
 ## load the data
 data("precipdata")

 ## Here the resolution of the data is in 'hours', we want to aggeregate the data up to 72 hours
 ## specify the aggregation durations

 durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)

 ## get the aggrageted data for each of the
 station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",  durations = durations)

}
