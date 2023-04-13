# quantiles of the IDF
q_egpd_idf = function(probs, par, dur) {
  sigma_d =  get_sigma_d(par = par, d = dur)
  quantiles = qextgp(p = probs, prob=NA, kappa=par[1], delta=NA, sigma= sigma_d, xi=par[3] , type=1)
  # if (xi_constant) {
  #   sigma_d =  ifelse(length(par)==4,par[2]/(dur)^par[4] ,par[2]/(dur+par[5])^par[4])
  #   #print(sigma_d)
  #   quantiles = qextgp(p = probs, prob=NA, kappa=par[1], delta=NA, sigma= sigma_d, xi=par[3] , type=1)
  # } else { # xi is not constant
  #   sigma_d =  ifelse(length(par)==3,par[2]/(dur)^par[3] ,par[2]/(dur+par[4])^par[3])
  #   #print(sigma_d)
  #   print(qqplot_egpd(data = station_data[seq(init_time_step,nrow(station_data),declustering_duration[id]),id], kappa = fit$par[1],
  #                     sigma = sigma_d, xi =  fitted_idf$xi_spline[id]))
  #   quantiles = qextgp(p = probs, prob=NA, kappa=par[1], delta=NA, sigma= sigma_d, xi=par[3] , type=1)
  # }
  names(quantiles) = "quantiles"
  return(quantiles)
}
#just a check
#function to get the non_exccedance probability for ecah duration d, depending on the return level
get_non_exc_probs_idf = function(x, y, npy, Tr_vec = c(2, 5, 10, 20, 50, 100), init_time_step=1){
  #x = station_data (df: col = data_d, row = observation)
  # y = vector of declustering durations
  # vector of empercial no of values per year that depends on the season for each data_d eg 90  for 24hours data in winter
  # Tr_vec: scalar or vector: Tr year or years intended

  data = x
  data = data[seq(init_time_step,length(data),y)]
  data = na.omit(data)
  n_postive_rain = length(data[data>0])
  n_year_obs = length(data)/npy
  probs = 1-1/(Tr_vec*n_postive_rain/n_year_obs)
  names(probs) = Tr_vec
  probs

}
#'  Quantile-Quantile Plot of EGPD model
#'
#'  The function plot qqplot from a fitted egpd model
#'
#' @param data vector of observation for a station
#' @param kappa  egpd 'kappa'  paramater
#' @param sigma  egpd 'sigma' paramater
#' @param xi    egpd shape parammter (xi)
#' @param ranking_method  a character string specifying how ties are treated. see \code{rank}  function in \code{R}. Defualts to "average"
#' @param title optional a character string for the title of the plot
#' @param subtitle optional a character string for the subtitle of the plot
#' @param confid logical, whether confidence bounds are to be plotted, defaults to \code{FALSE}
#' @param q_up  a numerical vector. For the upper confidence levels, only if \code{confid = T}
#' @param q_low  a numerical vector. For the lower confidence levels, only if \code{confid = T}
#' @details
#'   Note that  for confidence bounds to be plotted, \code{q_up} and \code{q_low} have to be separately computed by the user, eg via bootstrap
#'
#' @return a ggplot object
#'
#' @export
qqplot_egpd = function(data, kappa, sigma, xi, ranking_method = "average", title = "QQPLOT EGPD",
                       subtitle = "", confid = FALSE, q_up = NA, q_low = NA){

  data = na.omit(data)
  data = sort(data[data>0], decreasing = F)
  rank_i = rank(data, ties.method = ranking_method)
  q_emp  <- data.frame(d =  data, probs = rank_i/(length(data[data>0])+1))
  q_theo <- qextgp(p=q_emp$probs, prob=NA, kappa=kappa, delta=NA, sigma= sigma, xi=xi , type=1)
  data_gg = data.frame(x= q_emp$d, y = q_theo)
  max_q = max(data_gg$x, data_gg$y)
  qplot= ggplot(data_gg)+ geom_point( aes(x = x,y = y)) + geom_abline() + theme_void() +theme_bw() +
    scale_x_continuous(breaks = round(seq(0,30,5) * max_q/30, 1), limits = c(0, max_q)) +
    scale_y_continuous(breaks = round(seq(0,30,5) * max_q/30, 1), limits = c(0, max_q)) +
    labs(x="Empirical", y ="Theoretical", title = title, subtitle = subtitle)

  if (confid == T & !missing(q_up) & !missing(q_low)) {
    qplot = qplot +
      geom_line(aes(x = x, y=q.mle.U), col = "red") +
      geom_line(aes(x = x, y=q.mle.L), col = "red")
    #geom_ribbon(aes(x = x, ymin=q_low, ymax=q_up), alpha = 0.1)
  }
  return(qplot)
}

# tr plot
# function tp produce return level plot according to weibull equation
# Arguments: data : station series
#           kappa, sigma, xi: fitted egpd parameters
# returns: ggplot object of the return level plot

#'  Quantile-Quantile Plot of EGPD model
#'
#'  The function plot qqplot from a fitted egpd model
#'
#' @param data vector of observation for a station
#' @param kappa  egpd 'kappa'  paramater
#' @param sigma  egpd 'sigma' paramater
#' @param xi    egpd shape parammter (xi)
#' @param max_T  a scalar. The largest return level to be considered in the plot
#' @param npy number of observations per year. Eg if  \code{data} contains only daily data in summer, then \code{npy = 92}. If  for a year, then \code{npy = 365.25} etc description
#' @details
#'   to be done
#'
#' @return a ggplot object
#'
#' @export

trplot_egpd = function(data, kappa, sigma, xi, npy, max_T = 1000){
  data = na.omit(data)
  n_postive_rain = length(data[data>0])
  n_year_obs = length(data)/npy
  jj <- seq(-1, log10(max_T), by = 0.1)
  tr_data = data.frame(tr_years = 10^jj)
  tr_data$tr_quantiles <- NA
  for(i in 1:length(tr_data$tr_years)){
    tr_data$tr_quantiles[i] = qextgp(p=1-1/(tr_data$tr_years[i]*n_postive_rain/n_year_obs), prob=NA, kappa=kappa, delta=NA, sigma= sigma, xi=xi , type=1)
  }
  q_emp  <- data.frame(d = sort(data[data>0], decreasing = F) , probs = (1:length(data[data>0]))/(length(data[data>0])+1))
  #q_theo <- qextgp(p=q_emp$probs, prob=NA, kappa=kappa, delta=NA, sigma= sigma, xi=xi , type=1)
  data_emp = data.frame(x =  1/(1-q_emp$probs)*(n_year_obs/n_postive_rain), y= q_emp$d)
  tr_plot= ggplot()+ geom_point(data = data_emp, aes(x = x,y = y)) +
    geom_line(data = tr_data, aes(x = tr_years,y = tr_quantiles))+ scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                                                                 labels = scales::trans_format("log10", math_format(10^.x)))+   theme_void() +theme_bw() +
    labs(x="Return period (years)", y ="Return level (mm)", title = "Return level PLOT") + annotation_logticks(sides = "b")
  return(tr_plot)
}


#'  Plot IDF Curves
#'
#'  The function plot IDF curves from a fitted egpdIDF model
#'
#' @param station_data data_frame of observation for a station (nobs x aggregation duration)
#' @param kappa_fit  a numerical vector  of  egpd 'kappa'  parmater. Same length as \code{ncol(station_data)}
#' @param sigma_fit same as  \code{kappa_fit} but for egpd 'sigma'
#' @param xi_fit   same as  \code{kappa_fit} but for egpd shape parammter (xi)
#' @param durations  a (vector) of durations (same length as ncol of station_data)
#' @param declustering_duration a vector same length as \code{duration} .Whether the data should be temporally declustered, if yes the time step for each duration, or a vector of 1 otherwise
#' @param npy number of days per year. Eg if  \code{station_data} contains only data in summer, then \code{npy = 92}. If  for a year, then \code{npy = 365.25} etc
#' @param Tr_vec numerical vector of return periods fro which the curves are to be plotted.
#' @param init_time_step  a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd hour will be selected, and then a sequence is applied
#'  If a character, then a weighted normalized root mean square error will be used.
#'
#' @details
#'   to be added
#'
#' @return a ggplot object
#'
#' @export

plot_egpdidf_curves =function(station_data, kappa_fit, sigma_fit, xi_fit, durations, declustering_duration, npy, Tr_vec = c(2,5,10,20,50,100), init_time_step =1  ){

  #2.0 return level
  non_ex_probs = mapply(get_non_exc_probs_idf, x = station_data, y=  declustering_duration, npy = npy*24/declustering_duration,
                        init_time_step = init_time_step,MoreArgs = list(Tr_vec = Tr_vec) )

  # IDF_return_levels(station_data, declustering_duration , durations , Tr_vec, npy = npy, parameters = data.frame(kappa_fit, sigma_fit, xi_fit)) %>% data.frame %>%
  #    mutate(Tr = Tr_vec) %>%  tidyr::pivot_longer(cols = 1:length(durations),names_to = "dur", values_to = "tr_idf")
  #
  tr_levels = purrr::map_dfc(1:length(durations), function(i){
    qq= qextgp(p = non_ex_probs[,i], prob=NA, kappa =  kappa_fit[i],
               sigma= sigma_fit[i], xi=xi_fit[i], type=1)# * durations[i]
    return(setNames(data.frame(qq), durations[i]))
  }) %>% mutate(Tr = Tr_vec) %>%  tidyr::pivot_longer(cols = 1:length(durations),names_to = "dur", values_to = "tr_idf")
  emp_lev=emperical_levels_idf(station_data, declustering_duration ,non_ex_probs =non_ex_probs,  init_time_step = 1, durations ,Tr_year = Tr_vec,
                               npy = npy)  %>% data.frame() %>% setNames(durations)  %>% mutate(Tr = Tr_vec) %>%
    tidyr::pivot_longer(cols = 1:length(durations),names_to = "dur", values_to = "emp_lev")

  plot_data_idf = tr_levels %>% dplyr::mutate(emp_lev)

  p= ggplot()+
    geom_line(data = plot_data_idf, aes(x = as.numeric(dur),y = tr_idf, col = factor(Tr) ))+
    geom_point(data = plot_data_idf, aes(x = as.numeric(dur),y = emp_lev, col = factor(Tr))) +
    scale_x_continuous(trans = "log10", breaks = durations)+
    scale_y_continuous(trans = "log10")+theme_bw() +
    labs(x="Duration (hours)", y ="Return level (mm)", title = "IDF Curves", col = "Return Period") +
    #annotation_logticks(sides = "b") +
    #facet_wrap(.~factor(season, levels = seasons_name), scales = "free_y")  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  return(p)
}
##
#function to compute retirn levels
IDF_return_levels = function(station_data, declustering_duration, durations, Tr_year, parameters, npy, init_time_step=1){
  #npy: scalar in c(90,92,92,91)
  #Tr_year: sacalar or vector of return period eg 100 or  c(10,20,100)
  # parameters : dataframe of k,sigma, xi: each row corrspond to a particular duration

  # compute the non_exceedance prob for each TR and for each duration
  non_ex_probs = mapply(get_non_exc_probs_idf, x = station_data, y=  declustering_duration, npy = npy*24/declustering_duration, MoreArgs = list(Tr_vec = Tr_year, init_time_step=init_time_step) )


  if (is.vector(non_ex_probs))
    non_ex_probs = data.frame(t(non_ex_probs))

  return_levels_d = non_ex_probs
  # colnames(idf.array) = colnames(non_ex_probs)
  # rownames(idf.array)

  for (id in 1:length(durations)){
    return_levels_d[,id] = qextgp(p = non_ex_probs[,id], prob=NA, kappa = parameters[id,1], delta=NA, sigma= parameters[id,2], xi=parameters[id,3] , type=1)
  }

  return_levels_d
}
## same as above but returns emerical levels
emperical_levels_idf = function(station_data, declustering_duration, durations, Tr_year, npy = npy, non_ex_probs, init_time_step=1){

  #Tr_year: sacalar or vector of return period eg 100 or  c(10,20,100)
  # parameters : dataframe of k,sigma, xi: each row corrspond to a particular duration

  # compute the non_exceedance prob for each TR and for each duration
  if (missing(non_ex_probs)) {
    non_ex_probs = mapply(get_non_exc_probs_idf, x = station_data, y=  declustering_duration, npy = npy*24/declustering_duration, MoreArgs = list(Tr_vec = Tr_year, init_time_step=init_time_step) )

    if (is.vector(non_ex_probs))
      non_ex_probs = data.frame(t(non_ex_probs))
  }

  emperical_levels_d = non_ex_probs
  # colnames(idf.array) = colnames(non_ex_probs)
  # rownames(idf.array)

  for (id in 1:length(durations)){
    data = station_data[,id]
    data = data[seq(init_time_step,length(data),declustering_duration[id])]
    data = na.omit(data)
    emperical_levels_d[,id] = quantile(data[data>0], probs = non_ex_probs[,id], type = 6) #6 for the weibull

  }

  emperical_levels_d

}


#'  A wrapper function of \code{nrmse_idf}, computes NRMSE same for a range of durations
#'
#'  NRMSE stands for Normalized Root Mean Square Error
#'
#' @param station_data data_frame of observation for a station (nobs x aggregation duration)
#' @param kappa_fit  a numerical vector  of  egpd 'kappa'  parmater. Same length as \code{ncol(station_data)}
#' @param sigma_fit same as  \code{kappa_fit} but for egpd 'sigma'
#' @param xi_fit   same as  \code{kappa_fit} but for egpd shape parammter (xi)
#' @param declustering_duration a vector same length as \code{duration} .Whether the data should be temporally declustered, if yes the time step for each duration, or a vector of 1 otherwise
#' @param q optional. a number in percentage, or any character eg \code{nrsme_quantile = 90}, the normalized root mean square error will only be computed on excesses of the quantile. Quantile is over ALL  values (zeros included)
#'  If a character, then a weighted normalized root mean square error will be used.
#' @param init_time_step  a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd hour will be selected, and then a sequence is applied
#'
#' @details
#'   to be added
#'
#' @return a numerical vector of 'nrmse'
#'
#' @export
compute_nrsme = function(station_data, declustering_duration,  kappa_fit, sigma_fit, xi_fit,init_time_step=1,  q = 0){

  nrmse_d = sapply(1:ncol(station_data), function(id){
    res = NA
    try(res <- nrmse_idf(data = station_data[seq(init_time_step,nrow(station_data),declustering_duration[id]),id], kappa = kappa_fit[id],
                         sigma = sigma_fit[id], xi = xi_fit[id], q=q), silent = T)
    res
  })

  # for (id in 1:ncol(station_data)) {
  #   #d =  durations[id]
  #   #sigma_d[id] =   get_sigma_d(par = fit$par, d = d)
  #   nrmse_d[id] = nrmse_idf(data = station_data[seq(init_time_step,nrow(station_data),declustering_duration[id]),id], kappa = kappa_fit[id],
  #                           sigma = sigma_fit[id], xi = xi_fit[id], q=q)
  # }
  names(nrmse_d) = colnames(station_data)
  return(nrmse_d)
}
# function to predict from already fitted lienar IDF models. This is if the ftiited model is not availabe.
# The function fits a known relationship to the output of the original linear mmodel and then predicts from it
predict_from_linear_models =function(paramet, covar, model_type= 'log_log_MR', valid_covar){
  if(missing(valid_covar))
    valid_covar = covar
  init_breakponts = log(6)
  if (model_type == 'log_log_MR'){
    m4 = lm(log(paramet)~log(covar))
    y = log(paramet)
    x =log(covar)
    fit = m4
    #try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)) %>% suppressWarnings
    fitted_values =  predict.segmented(fit, newdata = data.frame(x = log(valid_covar))) %>%  exp

  } else if (model_type == 'linear_log_MR'){
    m4 = lm(paramet~log(covar))
    y = paramet
    x =log(covar)
    fit = m4
    #try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts))  %>% suppressWarnings
    fitted_values =  predict.segmented(fit, newdata = data.frame(x = log(valid_covar)))
    fitted_values[fitted_values<0] = 1e-6
  }
  names(fitted_values) = valid_covar
  return(fitted_values)
}
