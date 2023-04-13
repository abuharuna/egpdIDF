#' function to aggregate data into different intensities of longer durations
#'
#'  Given a precipitation data of particular duration, the function returns a dataframe containing the aggregate data
#'
#'
#' @param sample_data data_frame of at least two columns, one of which is  named \code{date} containg the dates, should be in correct \code{R} date format. The other column contains the timeseries for each station
#' @param durations  a numerical vector of aggregation durations. Must be in hours, if the resolution is hourly or in minutes. If in days, it should be in days, eg 1,2,3
#' @param st_code a character. Must be one of \code{colnames}  of \code{sample_data}.

#'
#' @details
#'   to be added
#'
#' @examples
#'  ## load the data
#'  data("precipdata")
#'
#'  ## Here the resolution of the data is in 'hours', we want to aggeregate the data up to 72 hours
#'  ## specify the aggregation durations
#'
#'  durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)
#'
#'  ## get the aggrageted data for each of the
#'  station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",  durations = durations)
#'
#' @return # A data frame: with \code{ncol = length(durations)}. Each containing the aggregated time series
#' @export
aggregate_data = function(sample_data, st_code, durations){

  list_df = lapply(st_code, function(x) data.frame(date=sample_data$date,RR=sample_data[, x]))
  names(list_df) = st_code
  ## check the time units of the data
  dtime <- ((list_df[[1]][2, "date"]- list_df[[1]][1,  "date"]))
  ds = durations
  if (attr(dtime, 'units') == "mins") #if true data is in minutes,  convert ds to minutes
    ds = ds*60

  rolled_data=roll(data = list_df, ds=ds)
  #df= data.frame(t(rolled_data))
  if (length(durations[durations < 1]) == 0) {
    if ((attr(dtime, 'units') == "days")) {
      colnames(rolled_data) = c(paste0("D_", c(durations*24), "h"))
    } else{
      colnames(rolled_data) = c(paste0("D_", c(durations), "h"))
    }

  } else {
    colnames(rolled_data) = c(paste0("D_", c(durations[durations < 1]*60), "m"),paste0("D_", c(durations[durations >= 1]), "h"))
  }

  #convert the intensities from "per 10mins" to "per hour", including  the initial 10mins
  if (attr(dtime, 'units') == "mins"){
    #if true data is in minutes,  convert ds to minutes
    station_data = rolled_data*6

  } else if (attr(dtime, 'units') == "days"){
    station_data = rolled_data/24
  }else{
    station_data = rolled_data
  }

  return(station_data)

}
