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
#
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
# function to plot the idf (MODIFY FOR TR INSTEAD OF P)
IDF.plot_egpd =function (durations, parameters, xi_fit, datadriven = F, probs = c(0.99, 0.999, 0.9999), title= "IDF_egpd", cols = 4:2, ylim,xlabs,ylabs, ...) {
  if (datadriven & !inherits(parameters, "data.frame"))
    stop("You have specfied  'datadriven = T', the argument 'parameters', should be a dataframe of kappa, sigma and xi")
  if (is.vector(probs) ) {
    if (length(cols) != length(probs)) cols <- 1:nrow(probs)
  } else {
    if (length(cols) != nrow(probs)) cols <- 1:nrow(probs)
  }

  # to store the levels
  idf.array = matrix(data = NA, nrow(probs), length(durations))
  if (!datadriven & any(names(parameters) == "xi0")) { #constant xi
    # qs <- lapply(durations, q_egpd_idf, probs = probs,  par=parameters)
    # idf.array <- matrix(unlist(qs), length(probs), length(durations))

    #idf.array = matrix(data = NA, length(probs), length(durations))
    par = parameters
    for (id in 1:length(durations)){
      dur = durations[id]
      sigma_d =  get_sigma_d(par = par, d = dur)
      #sigma_d =  ifelse(length(par)==3,par[2]/(dur)^par[3] ,par[2]/(dur+par[4])^par[3])
      idf.array[,id ] = qextgp(p = probs[,id], prob=NA, kappa =par[1], delta=NA, sigma= sigma_d, xi=par['xi0'] , type=1)
    }

  } else if(!datadriven) { # xi function of duration
    #idf.array = matrix(data = NA, length(probs), length(durations))
    par = parameters
    for (id in 1:length(durations)){
      dur = durations[id]
      sigma_d =  get_sigma_d(par = par, d = dur)
      #sigma_d =  ifelse(length(par)==3,par[2]/(dur)^par[3] ,par[2]/(dur+par[4])^par[3])
      idf.array[,id ] = qextgp(p = probs[,id], prob=NA, kappa =par[1], delta=NA, sigma= sigma_d, xi=xi_fit[id] , type=1)
    }
  } else if (datadriven){
    #idf.array = matrix(data = NA, length(probs), length(durations))
    for (id in 1:length(durations)){
      idf.array[,id] = qextgp(p = probs[,id], prob=NA, kappa =parameters[id,1], delta=NA, sigma= parameters[id,2], xi=parameters[id,3] , type=1)
    }
  }


  if (missing(ylim))
    ylim = c(min(idf.array), max(idf.array))

  if (missing(xlabs))
    xlabs = "Duration"

  if (missing(ylabs))
    ylabs =  "Intensity [mm/hr]"

  plot(NA, xlim = c(min(durations), max(durations)), ylim = ylim,  xlab = xlabs, ylab = ylabs,  log = "xy", main = title, ...)
  for (i in 1:nrow(probs)) {
    lines(durations, idf.array[i, ], col = cols[i],  lwd = 2)
  }
  legend(x = "topright", ncol = 1 , title = "T-years", legend = rev(rownames(probs)), fill = rev(cols),  bty = 'n', cex = 0.8 )

}

##
#function to compute retirn levels
IDF_return_levels = function(station_data, declustering_duration, durations, Tr_year, parameters, npy, init_time_step=1){
  #npy: scalar in c(90,92,92,91)
  #Tr_year: sacalar or vector of return period eg 100 or  c(10,20,100)
  # parameters : dataframe of k,sigma, xi: each row corrspond to a particular duration

  # compute the non_exceedance prob for each TR and for each duration
  non_ex_probs = mapply(get_non_exc_probs_idf, x = station_data, y=  declustering_duration, npy = npy*24/durations, MoreArgs = list(Tr_vec = Tr_year, init_time_step=init_time_step) )


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
    non_ex_probs = mapply(get_non_exc_probs_idf, x = station_data, y=  declustering_duration, npy = npy*24/durations, MoreArgs = list(Tr_vec = Tr_year, init_time_step=init_time_step) )

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
#eg call
#IDF_return_levels(station_data, declustering_duration , durations , c(10,100), npy = 90, parameters = data.frame(kappa_fit, sigma_fit, xi_fit))
#emperical_levels_idf(station_data, declustering_duration , durations , c(10,100), npy = 90, parameters = data.frame(kappa_fit, sigma_fit, xi_fit))




#wrapper function to nrmse_idf, computes same for a reange od durations
#q is in % eg 0% for all postive data, and 95% for exceedances of a 95% quantile. Quantile is over ALL  values (zeros included)
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
