


#'  Function to estimates the initial parameters for the IDF_egpd models
#'
#'  each duration is fitted separately, then lm is fitted to estimate the parameters of the IDF (sigma0 as the intecept, etha as the -ve of the slope)
#'  kappa and xi as the median of the d_ estimates
#' @param station_data data_frame of observation for a station (nobs x aggregation duration)
#' @param durations  a (vector) of durations (same length as ncol of station_data)
#' @param declustering_duration a vector same length as \code{duration} .Whether the data should be temporally declustered, if yes the time step for each duration, or a vector of 1 otherwise
#' @param init_time_step  a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd hour will be selected, and then a sequence is applied
#' @param fitting_method eiither \code{"mle"} for maximum likelihood or \code{"pwm"} for probability weighted momoments. Defaults to \code{"mle"}
#' @param auto_fit logical, whether an automatic fit is to be done to find a censoring threshold that minimizes the normalized root mean square error; defaults to \code{TRUE}
#' @param nrmse_tol a scalar, defaults to \code{0.1}. defaines a relatively good egpd fit. If higher, and \code{auto_fit=TRUE}, automatic fit is done
#' @param simple_scaling logical, defaults to \code{TRUE}, defines the initial values values that would be returned
#' @param use_r_optim logical, defaults  to \code{FALSE}, if yes, \code{optim} will be used to find the best lower censoring values, otherwise, descrete values will be tested
#' @param nrsme_quantile a number in percentage, or any character eg \code{nrsme_quantile = 90}, the normalized root mean square error will only be computed on excesses of the quaneile.  If a charanter, then a weighted normalized root mean square error will be used.
#'
#' @details
#'   to be added
#'
#' @return # A list:
#'    init: the initial parameters to use for fitting a egpd IDF model
#'    fits: egdp parameters fit for each duration seperately
#'    regres_summary: summary of the lm fit of the sigma vs duration
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
#' \dontrun{
#'  initial_params = egpd_idf_init(station_data = station_data,
#'                 durations = durations, fitting_method = "mle",
#'                 declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72), auto_fit = FALSE)
#'  ## check the fitted egpd parameters for each duration
#'  initial_params$fits$kappa_param
#'  initial_params$fits$scale_param
#'  initial_params$fits$shape_param
#'
#'  ## check the quality of the fit
#'  initial_params$fits$nrsme
#'
#'  ## for a good fit, 'nrsme" should be small.
#'
#'  ##  its always good to  use left censoring with 'mle' fit. Lets try, and check the 'nrmse' again
#'  ##  we set "auto_fit=T", "nrmse_tol=0.1"  "use_r_optim=T" , "nrsme_quantile = 0".
#'  ## Check the arguments for their meaning
#'  initial_params = egpd_idf_init(station_data = station_data,
#'            durations = durations, fitting_method = "mle",
#'             declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
#'               auto_fit = TRUE, nrmse_tol = 0.1,use_r_optim = TRUE, nrsme_quantile = 0)
#'  ## check the quality of the fit
#'  initial_params$fits$nrsme
#'  ## check the parameters
#'  initial_params$fits$kappa_param
#'  initial_params$fits$scale_param
#'  initial_params$fits$shape_param
#' }
#' @export
egpd_idf_init <- function(station_data, durations, declustering_duration, init_time_step = 1,
                          fitting_method = "mle", auto_fit = T, nrmse_tol = 0.1,simple_scaling = T, use_r_optim = F, nrsme_quantile = 0){

  # fit the data for each duration seperately
  if (use_r_optim) { # r optim function used to find the best lower_censoring; slower but find the best threshold
    station_fit <- local_fit_IDF_h_par(sample = station_data, fitting_method= fitting_method, auto_fit = auto_fit, nrmse_tol = nrmse_tol,
                                       low_th = c(seq(0,1.2,0.1)), durations = durations, declustering_duration = declustering_duration,
                                       init_time_step=init_time_step, q = nrsme_quantile)

  } else{ # here discrete values of lower threshold are tested, its possible to miss the best value
    station_fit <- local_fit_IDF_h(sample = station_data, fitting_method= fitting_method, auto_fit = auto_fit, nrmse_tol = nrmse_tol,
                                   low_th = c(seq(0,1.2,0.1)), durations = durations, declustering_duration = declustering_duration,
                                   init_time_step=init_time_step, q = nrsme_quantile)

  }


  #get initial params for scaling model
  obj = prepare_idf_obj(fit_obj = station_fit, simple_scaling = simple_scaling, durations = durations)


  return(append(obj, list('fits' = station_fit)))
  #return(list("init" = init, "fits" = station_fit, "regres_summary" = summary(lmsig)))
}

#function to get initila paramters for IDF model
prepare_idf_obj <-function(fit_obj, simple_scaling = T, durations){

  # fit_obj: objected returned by local_fit of data of diffrent duration
  #extract the fitted parpams
  sigma_vec = fit_obj$scale_param
  kappa_vec = fit_obj$kappa_param
  xi_vec = fit_obj$shape_param

  # estimaete the params of the IDF
  # by definiton, only sigma of EGPD shows scaling
  lmsig <- lm(log(sigma_vec)~log(durations))
  sigma0 = exp(coefficients(lmsig)[1]) # intercept of the fit
  eta0 =  -(coefficients(lmsig)[2]) # -ve slope of the fit

  # kappa and xi should be uniform (indepernedent of duration)
  kappa0 = median(kappa_vec)
  xi0 = median(xi_vec)

  #initial theta
  theta0 = 0

  init = c(kappa0, sigma0, xi0, eta0)
  names(init) = c("kappa0", "sigma0", "xi0", "eta")

  #if not simple scaling, then we need to include theta
  if (!simple_scaling) {
    init[5] = theta0
    names(init)[5] = "theta"
  }
  list("init" = init,  "regres_summary" = summary(lmsig))
}


### For hourly resolution
local_fit_IDF_h  <- function(sample, fitting_method= "mle", auto_fit = T, nrmse_tol = 0.1, low_th = c(seq(0,1,0.1)), durations,
                             declustering_duration, init_time_step=1, q = 0){
  starting_c = 0#ifelse(fitting_method== "mle", 0.5,0)
  scale_param <- c()
  shape_param <- c()
  kappa_param <- c()
  rmse <- c()
  lower_c <- rep(starting_c, ncol(sample))/durations
  rounding_c <- rep(0, ncol(sample))/durations
  cens_thresholds_init = expand.grid(low_th = low_th, rounding= c(0))
  for (i in 1:ncol(sample)) {
    data <- sample[seq(init_time_step,nrow(sample),declustering_duration[i]),i]  #temporal declustering
    data = na.omit(data)
    data = data[data<2000] #omit  Huge outliers
    inits = fevd(data, method = "MLE", type="GP", threshold=0, optim.args=list(method="L-BFGS-B", lower =c(0.1, 0.0001), upper = c(40, 0.3))  )$results$par

    extgp <- fit.extgp(data = data[data>0], model=1, method = fitting_method, init =  c(0.9, inits) ,
                       censoring=c(starting_c/durations[i],Inf), rounded = rounding_c[i], confint = F,  plots = F, R = 1)
    kappa_param[i] <- ifelse(fitting_method=="mle", extgp$fit$mle[1], extgp$fit$pwm[1] )
    scale_param[i] <- ifelse(fitting_method=="mle", extgp$fit$mle[2], extgp$fit$pwm[2] )
    shape_param[i] <- ifelse(fitting_method=="mle", extgp$fit$mle[3], extgp$fit$pwm[3] )
    # nrmse computation
    #**** normalized rmse to track quality of qqplot *******#
    #sorted_data = sort(data[data>0], decreasing = F)
    #empirical <- data.frame(d = sorted_data, probs = (rank(sorted_data, ties.method = "average"))/(length(data[data>0])+1))
    #mle = qextgp(p=empirical$probs, prob=NA, kappa= kappa_param[i], delta=NA, sigma=scale_param[i], xi= shape_param[i], type=1)
    #top_5 <- ceiling(length(empirical$d)*.95):length(empirical$d)
    #rmse[i] = sqrt(mean((mle - empirical$d)^2))/mean(data[data>0], na.rm=T)
    rmse[i] = nrmse_idf(data ,kappa= kappa_param[i], sigma=scale_param[i], xi= shape_param[i], q = q)
    cens_thresholds = cens_thresholds_init/durations[i]
    #* if autofit, then we vary the censoring thresholds and compute nrmse each time. finally we retain the params with the best fit (least nrmse)
    if (auto_fit & rmse[i] >  nrmse_tol) {
      #message("auto_fit, searching for best censoring thresholds")
      new_rmse = c()
      new_params = data.frame(array(NA, dim = c(nrow(cens_thresholds), 3)))
      colnames(new_params) = c("kappa", "sigma", "xi")
      for (r in 1:nrow(cens_thresholds)) {
        try_fit = fit.extgp(data = data[data>0], model=1, method = fitting_method, init =  c(0.9, inits) ,
                            censoring=c(cens_thresholds[r,1],Inf), rounded = cens_thresholds[r,2], confint = F,  plots = F)
        if (fitting_method=="mle")  new_params[r,] = try_fit$fit$mle
        if (fitting_method=="pwm")  new_params[r,] = try_fit$fit$pwm
        #new_params[r,] = ifelse(, try_fit$fit$mle, try_fit$fit$pwm )
        # rmse computation
        #mle = qextgp(p=empirical$probs, prob=NA, kappa= new_params$kappa[r], delta=NA, sigma=new_params$sigma[r], xi= new_params$xi[r], type=1)
        #new_rmse[r] = sqrt(mean((mle - empirical$d)^2))/mean(data[data>0], na.rm=T)
        new_rmse[r]  = nrmse_idf(data ,kappa= new_params$kappa[r],  sigma=new_params$sigma[r], xi= new_params$xi[r], q = q)
      }
      best_set = which.min(new_rmse)
      kappa_param[i] <- new_params$kappa[best_set]
      scale_param[i] <- new_params$sigma[best_set]
      shape_param[i] <- new_params$xi[best_set]
      rmse[i] <- new_rmse[best_set]
      rounding_c[i] <- cens_thresholds$rounding[best_set]
      lower_c[i] <- cens_thresholds$low_th[best_set]
    }
    # n_year_obs <- length(data)/npy
    # returns100[i] <- qextgp(p=(1-1/(rp*length(data[data>0])/n_year_obs)), prob=NA, kappa= kappa_param[i], delta=NA, sigma=scale_param[i], xi= shape_param[i], type=1)
    #message(paste0("station ", i, " out of ", length(sample)))
  }
  names(scale_param) = names(shape_param) = names(kappa_param) = colnames(sample)
  return(list("sample" = sample, "scale_param" = scale_param, "shape_param" = shape_param, "kappa_param"= kappa_param,
              "nrsme" = rmse, "lower_threshold"= lower_c, "rounding" = rounding_c))
}
# Function linked to the automatic lower censoring threshold determination in the EGPD fit
#this fucntion computes the NRMSE corresponding to a particular choice of threshold.
# it is minimized by the optim fucntion called in the **local_fit_IDF_h_auto** function above
auto_egpd_fit = function(data, fitting_method , inits, lower_cens, q = 0){
  # fit the parameters
  data = na.omit(data)
  par.fit = fit.extgp(data = data[data>0], model=1, method = fitting_method, init =  c(0.9, inits) ,
                      censoring=c(lower_cens,Inf), rounded = 0, confint = F,  plots = F, R = 1)$fit[[1]]
  # rmse computation

  # sorted_data = sort(data[data>0], decreasing = F)
  # empirical <- data.frame(d = sorted_data, probs = (rank(sorted_data, ties.method = "average"))/(length(data[data>0])+1))
  #
  # mle = qextgp(p=empirical$probs, kappa = par.fit[1], sigma = par.fit[2], xi = par.fit[3], type=1)
  # rmse=sqrt(mean((mle - empirical$d)^2))/mean(data[data>0], na.rm=T)

  rmse = nrmse_idf(data ,kappa = par.fit[1], sigma = par.fit[2], xi = par.fit[3], q = q)
  return(rmse)
}
# fitti

#local fit of data with intensity per hour
# for each packge used, r optim used to fine the best lower cenosirng threholds
local_fit_IDF_h_par  <- function(sample, fitting_method= "mle", auto_fit = T, nrmse_tol = 0.05,  low_th = c(seq(0,1,0.1)),
                                 durations, declustering_duration, init_time_step=1, q = 0){
  #rmse = c()
  starting_c =0# ifelse(fitting_method== "mle", 0.5,0)
  #cens_thresholds = expand.grid(low_th = c(1,2,5), rounding= c(0.1, 1, 2, 2.5))
  n.cores <- min(5, length(durations))
  cl <- makeCluster(n.cores)
  clusterExport(cl, c('auto_egpd_fit','nrmse_idf') )
  registerDoParallel(cl)
  # start= Sys.time()
  fits= foreach(i = 1:ncol(sample), .packages = c("extRemes", "mev"), .errorhandling = "remove") %dopar% {
    data <- sample[seq(init_time_step,nrow(sample),declustering_duration[i]),i]
    data = na.omit(data)
    inits = fevd(data, method = "MLE", type="GP", threshold=0, optim.args=list(method="L-BFGS-B", lower =c(0.1, 0.0001), upper = c(100, 0.1))  )$results$par
    lower_c = starting_c/durations[i]
    extgp <- fit.extgp(data = data[data>0], model=1, method = fitting_method, init =  c(0.9, inits) ,
                       censoring=c(lower_c,Inf), rounded = 0, confint = F,  plots = F, R = 1)
    kappa <-  extgp[[1]][[1]][1]
    scale <-  extgp[[1]][[1]][2]
    shape <-  extgp[[1]][[1]][3]
    # nrmse computation
    #*normalized rmse to track quality of qqplot *******#
    #sorted_data = sort(data[data>0], decreasing = F)
    #empirical <- data.frame(d = sorted_data, probs = (rank(sorted_data, ties.method = "average"))/(length(data[data>0])+1))
    #mle = qextgp(p=empirical$probs, prob=NA, kappa= kappa, delta=NA, sigma=scale, xi= shape, type=1)
    #top_5 <- ceiling(length(empirical$d)*.95):length(empirical$d)
    #rmse = sqrt(mean((mle - empirical$d)^2))/mean(data[data>0], na.rm=T)
    rmse = nrmse_idf(data ,kappa= kappa, sigma=scale, xi= shape, q = q)
    #* if autofit, then we vary the censoring thresholds and compute nrmse each time. finally we retain the params with the best fit (least nrmse)
    if (auto_fit & rmse>  nrmse_tol) {
      result = optim(par = 1, fn = auto_egpd_fit, gr = NULL, data = data, fitting_method = fitting_method, q = q,
                     inits = inits,  method = "Brent", lower = min(data[data>0]),
                     upper = quantile(data[data>0], 0.9)) #ifelse(quantile(data[data>0], 0.9)>15,15,quantile(data[data>0], 0.9) )
      par.fit = fit.extgp(data = data[data>0], model=1, method = fitting_method, init =  c(0.9, inits) ,
                          censoring=c(result[[1]],Inf), rounded = 0, confint = F,  plots = F, R = 1)$fit[[1]]
      kappa <-  par.fit[1]
      scale <-  par.fit[2]
      shape <-  par.fit[3]
      rmse  <- result[[2]]
      #rounding_c <- 0
      lower_c <- result[[1]]
    }
    # if(rmse_o < rmse) { #check of the optimwed rmse is better than the initial  nrmse
    #   result_f = list("scale_param" = scale_o, "shape_param" = shape_o, "lower_threshold" = lower_c_o, "nrsme" = rmse_o, "kappa_param"= kappa_o)
    # } else{
    #   result_f = list("scale_param" = scale, "shape_param" = shape, "lower_threshold" = lower_c, "nrsme" = rmse, "kappa_param"= kappa)
    # }
    # result_f
    #message(paste0("station ", i, " out of ", length(sample)  ))
    list("scale_param" = scale, "shape_param" = shape, "lower_threshold" = lower_c, "nrsme" = rmse, "kappa_param"= kappa)
  }
  #Sys.time() - start
  stopCluster(cl)
  scale_param <- sapply(fits, function(y) y$scale_param)
  shape_param <- sapply(fits, function(y) y$shape_param)
  kappa_param <- sapply(fits, function(y) y$kappa_param)
  lower_threshold <- sapply(fits, function(y) y$lower_threshold)
  rmse <- sapply(fits, function(y) y$nrsme)

  names(scale_param) = names(shape_param) = names(kappa_param) = colnames(sample)
  return(list("sample" = sample, "scale_param" = scale_param, "shape_param" = shape_param,"lower_threshold" = lower_threshold,
              "nrsme" = rmse,  "kappa_param"= kappa_param))
}



#' function for fitting of data-driven IDF models
#'
#'  a 10 parameter IDF model in which the form of the EGPD parameter dependence is emperically
#'  determined from the data it self  For details, see Haruna et al 2023
#'
#'
#' @param station_data data_frame of observation for a station (nobs x aggregation duration)
#' @param durations  a (vector) of durations (same length as ncol of station_data)
#' @param declustering_duration a vector same length as \code{duration} .Whether the data should be temporally declustered, if yes the time step for each duration, or a vector of 1 otherwise
#' @param initial_params an object returned by \code{egpd_idf_init}.
#' @param censored a scalar or vector of \code{length(durations)}, for the left censoring to be applied to data of each duration. If a scalar is given, it will be divided by \code{durations}
#' @param init_time_step  a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd hour will be selected, and then a sequence is applied
#' @param fitting_method either \code{"mle"} for maximum likelihood or \code{"pwm"} for probability weighted momoments. Defaults to \code{"mle"}
#' @param use_mle_init logical, defaults  to \code{FALSE}, if yes, an iterative pairwise likelohoof fitting is done. See ...
#' @param optim_algo the \code{optim} algtorthm to use. defaults to \code{"Nelder-Mead" }
#'
#' @details
#'   to be added
#'
#' @return # A list:
#'
#'    fits: egdp parameters fit for each duration seperately
#'    optim details, to check convergence
#' @examples
#'  ## load the data
#'  data("precipdata")
#'
#'  ## Here the resolution of the data is 'hours', we want to aggeregate the data up to 72 hours
#'  ## specify the aggregation durations
#'
#'  durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)
#'
#'  ## get the aggrageted data for each of the
#' station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",
#'  durations = durations)
#' \dontrun{
#'  ## get initial values
#'
#'  initial_params = egpd_idf_init(station_data = station_data,
#'              durations = durations, fitting_method = "mle",
#'                declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
#'                auto_fit = T, nrmse_tol = 0.1,use_r_optim = T, nrsme_quantile = 0)
#'  ## fit the data driven IDF
#'  fitted_idf = fit_egpd_idf_data_driven(station_data = station_data, durations = durations,
#'              declustering_duration = c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72),
#'              fitting_method = 'mle',  initial_params = initial_params,
#'               optim_algo = "BFGS")
#'
#'  #check 'optim' params, for convergencem etc
#'  fitted_idf$fitted_params
#'
#'  #parameters of the IDF
#'  fitted_idf$fitted_params$par
#'
#'  # fitted egpd parameters for the given durations
#'  kappa_fit =  fitted_idf$kappa_param
#'  sigma_fit = fitted_idf$scale_param
#'  xi_fit = fitted_idf$shape_param
#'
#'  #compute nrmse to check quality of fit
#'  nrmse_d = compute_nrsme(station_data, c(1,2,3,6,10,12, 16,18, 24, 48, 72),
#'   kappa_fit, sigma_fit, xi_fit, init_time_step = 1, q = 0)
#'  nrmse_d
#'
#' # Plot the IDF curves
#' plot_egpdidf_curves(station_data = station_data,  kappa_fit = kappa_fit,
#' sigma_fit = sigma_fit, xi_fit = xi_fit, durations,
#'  declustering_duration=c(1,2,3,6,10,12, 16,18, 24, 48, 72), npy = 92, init_time_step=1 )
#'  }

#' @export
fit_egpd_idf_data_driven <- function(station_data, durations, declustering_duration, initial_params,
                                         censored,   fitting_method = "mle", init_time_step=1,
                                        use_mle_init =F,  optim_algo = "Nelder-Mead"){




  if (missing(initial_params)) {stop("provide initial paramters")} else{
    init =  initial_params$init
    #init = c(init[1],0,log(init[2]),0, init[3], 0)
    station_fit = initial_params$fits
  }

  if (missing(censored)){
    censored = initial_params$fits$lower_threshold
  } else if (length(censored)==1){
    censored = censored/durations
  }


  #xxx = fit_segmented(paramet = initial_params$fits$kappa_param, covar = durations)
  init_kappa = find_init_linear_models(paramet = initial_params$fits$kappa_param, covar = durations, model_type = 'log_log_MR', param_name = 'kappa')
  init_sigma = find_init_linear_models(paramet = initial_params$fits$scale_param, covar = durations, model_type = 'log_log_MR', param_name = 'sigma')
  #init_kappa = lm(formula = initial_params$fits$kappa_param ~ durations)$coefficients
  init_xi = lm((initial_params$fits$shape_param)~log(durations))$coefficients
  names(init_xi) = c('xi0', 'xi_eta1')

  scaling_breaks =  c(init_kappa[4], init_sigma[4])
  init = c(init_kappa[-4], init_sigma[-4], init_xi)

  #preparae tha data into a suitable format
  n_cens = map_dbl(seq_along(declustering_duration), function(d){

    data_dec =  decluster_set_NA(x = station_data[[d]], init_time_step = init_time_step, step = declustering_duration[d]) %>% na.omit
    # censoring_value = init_param %>%  filter(duration == durations[d], area == grid_area[a]) %>%  pull(lower_c)
    length(data_dec[data_dec < censored[d]])
  })


  station_data_censored = map_dfc(seq_along(declustering_duration), function(d){
    data_dec =  decluster_set_NA(x = station_data[[d]], init_time_step = init_time_step, step = declustering_duration[d])

    # seclect obs that are >= censoring value for duration d
    #censoring_value = init_param %>%  filter(duration == durations[d], area == grid_area[a]) %>%  pull(lower_c)
    data_dec[data_dec < censored[d]] =  NA
    data_dec = data.frame( data_dec ) %>%  setNames(durations[d])
    data_dec
  }) %>% #pivot longer
    tidyr::pivot_longer(cols = 1:length(durations),names_to = "duration", values_to = "precip",
                        names_transform = list(duration = as.numeric) ) %>%  arrange(duration)

  station_data_censored = tidyr::drop_na(station_data_censored)

  censored_data = station_data_censored %>%  group_by(duration) %>% group_split()
  #init = c(init_kappa, log(xxx['sigma01']), -xxx['eta1'], log(xxx['sigma02']), -xxx['eta2'] , init_xi, xxx['k'])
  #names(init) = c("kappa0", "kappa_slope", "sigma01",  "eta1" ,  "sigma02",  "eta2", "xi0", "xi0_slope",  "k" )
  # function for the likelihhod (with Xi dependent on D) to be mqximized taking into acoount lower censoring
  LK.EGPD.idf_GLM_new <- function(init, free_params, scaling_breaks ,censored_data,  n_cens, durations, censored){
    #### INPUT
    # init (vector(5)): kappa; sigma (scale), xi (shape), eta, theta (IDF of koutsoyiannis, 1998)
    # station_data: data_frame of observation for a station (nobs x aggreation duration)
    # durations: (vector) aggeragtion durations (same lenght as ncol of station_data)
    #censored: (scalar) censoring threshold (lower) to be applied to all durations ()
    # xi_fit : predicted xi from the smoothing spline regression

    #### OUTPUT
    # Likelihood value
    # init[1:2] = kappalink$linkfun(init[1:2])
    # init[3:6] = sigmalink$linkfun(init[3:6])
    # init[7:8] = xilink$linkfun(init[7:8])
    free_params[names(init)] = init

    kappa_init = free_params[which(substr(names(free_params), 1, 2) == "ka")]
    sigma_init = free_params[which(substr(names(free_params), 1, 2) == "si")]
    xi_init = free_params[which(substr(names(free_params), 1, 2) == "xi")]

    #add the scaling braeks
    kappa_init = c(kappa_init, scaling_breaks["kappa_k"])
    sigma_init = c(sigma_init, scaling_breaks["sigma_k"])

    #add thz second intercept
    kappa_init['kappa02'] = get_second_intercept(kappa_init)
    sigma_init['sigma02'] = get_second_intercept(sigma_init)


    vec_kappa_d = map_dbl(.x = durations, get_vector_of_param2, init_param = kappa_init)
    vec_sigma_d = map_dbl(.x = durations, get_vector_of_param2, init_param = sigma_init)
    vec_xi_d = map_dbl(.x = durations, ~  c(1, log(.x)) %*% xi_init)

    vec_kappa_d_full = sapply(seq_along(durations), function(i) rep(vec_kappa_d[i], length(censored_data[[i]]$precip)) )%>% unlist
    vec_sigma_d_full = sapply(seq_along(durations), function(i) rep(vec_sigma_d[i], length(censored_data[[i]]$precip)) )%>% unlist
    vec_xi_d_full = sapply(seq_along(durations), function(i) rep(vec_xi_d[i], length(censored_data[[i]]$precip)) )%>% unlist
    #kappa_vec = rep(init[1], nobs)
    #xi_vec =  rep(init[3], nobs)
    #vec_xi_d[vec_xi_d<1e-6] = 1e-6
    # if (all(vec_xi_d < 0) |  any(vec_xi_d > 0.5) ) {
    #   return(1e100)
    # }else{
    #   vec_xi_d[vec_xi_d<1e-6] = 1e-6
    # }
    if (all(vec_xi_d < 0) |  any(vec_xi_d > 0.5) ) {
      return(1e100)
    }else{
      vec_xi_d[vec_xi_d<1e-6] = 1e-6
    }

    # do a likelihood fit (with left censoring depending on the censored value set to zero or otherwise)
    # likelhood

    if(all(c(vec_kappa_d, vec_sigma_d)> 0)){
      # for left censored data
      is_censored = n_cens > 0
      contrib.cens1 <- ifelse(any(censored>0), sum(n_cens[is_censored] * log(p_egpd(x = censored[is_censored], kappa = vec_kappa_d[is_censored],
                                                                                    sigma_d  = vec_sigma_d[is_censored], xi = vec_xi_d[is_censored]))), 0)
      #for uncensored data
      # contrib.not.cens <- map_dbl(seq_along(vec_sigma_d), function(i){
      #   sum(d_egpd(x = censored_data[[i]]$precip, kappa =rep(vec_kappa_d[i], length(censored_data[[i]]$precip)),
      #              sigma_d = rep(vec_sigma_d[i], length(censored_data[[i]]$precip)),
      #              xi = rep(vec_xi_d[i], length(censored_data[[i]]$precip))), na.rm = T)
      #   #sum(dextgp(x = censored_data[[i]]$precip, kappa = vec_kappa_d[i], sigma = vec_sigma_d[i], xi = vec_xi_d[i], log = T), na.rm = T)
      # }) %>%  sum
      #sum(dextgp(x = censored_data[[i]]$precip, kappa = vec_kappa_d[i], sigma = vec_sigma_d[i], xi = vec_xi_d[i], log = T), na.rm = T)
      contrib.not.cens <- sum(d_egpd(x = sapply(seq_along(durations), function(i) censored_data[[i]]$precip) %>% unlist,
                                     kappa = vec_kappa_d_full, sigma_d = vec_sigma_d_full, xi = vec_xi_d_full))
      #print(-(contrib.cens1+ contrib.not.cens))
      #sum the two likelihhods
      return(-(contrib.cens1+ contrib.not.cens))
    }else{
      return(1e100)
    }
  }
  # LK.EGPD.idf_GLM(init =init, scaling_breaks = scaling_breaks, station_data=station_data, durations = durations, censored= censored,
  #                 declustering_duration = declustering_duration, init_time_step=init_time_step )

  # LK.EGPD.idf_GLM_new(init = init, free_params = init, scaling_breaks = scaling_breaks, censored_data = censored_data, n_cens = n_cens, durations = durations, censored = censored )
  #
  # # # optimize the function above
  # #message("IDF fitting")
  # par.optim = optim(par=init,fn=LK.EGPD.idf_GLM,gr=NULL,free_params = init,scaling_breaks = scaling_breaks, censored_data = censored_data,
  #                   n_cens = n_cens, durations = durations, censored = censored, control = list(maxit = 3000), hessian = FALSE,method="Nelder-Mead")


  if (use_mle_init) {
    i = 1
    tol = 10000
    likl = LK.EGPD.idf_GLM_new(init =init, free_params = init,  scaling_breaks = scaling_breaks, censored_data = censored_data,
                               n_cens = n_cens, durations = durations, censored = censored)
    init_temp = init
    while (tol > 0.05) {
      message(i)
      for (pp in c("ka", "si", "xi")) {
        pname =  which(substr(names(init), 1, 2) == pp)
        par.optim = optim(par=init_temp[pname], fn=LK.EGPD.idf_GLM_new, gr=NULL, free_params = init_temp, scaling_breaks = scaling_breaks, censored_data = censored_data,
                          n_cens = n_cens, durations = durations, censored = censored, control = list(maxit = 3000), hessian = FALSE,method=optim_algo)

        init_temp[pname] = par.optim$par
        print(c(par.optim$value, par.optim$convergence))
      }
      tol = likl -  par.optim$value
      print(tol)
      likl = par.optim$value
      i =i +1
    }
    init_ml = init_temp



    par.optim =  optim(par=init_ml, fn=LK.EGPD.idf_GLM_new, gr=NULL, free_params = init_ml, scaling_breaks = scaling_breaks, censored_data = censored_data,
                       n_cens = n_cens, durations = durations, censored = censored, control = list(maxit = 3000), hessian = FALSE,method="BFGS")
  } else {
    par.optim =  optim(par=init,fn=LK.EGPD.idf_GLM_new,gr=NULL,free_params = init, scaling_breaks = scaling_breaks,  censored_data = censored_data,
                       n_cens = n_cens, durations = durations, censored = censored, control = list(maxit = 3000), hessian = FALSE,method=optim_algo)
  }
  par = par.optim$par
  #prdictions
  kappa_par = par[which(substr(names(par), 1, 2) == "ka")]
  sigma_par = par[which(substr(names(par), 1, 2) == "si")]
  xi_par = par[which(substr(names(par), 1, 2) == "xi")]

  #add the scaling breaks
  kappa_par = c(kappa_par, scaling_breaks["kappa_k"])
  sigma_par = c(sigma_par, scaling_breaks["sigma_k"])

  #add the second intercept
  kappa_par['kappa02'] = get_second_intercept(kappa_par)
  sigma_par['sigma02'] = get_second_intercept(sigma_par)
  n_covariates = length(durations)
  #kappa
  #kappa_param = (cbind(rep(1, n_covariates), durations) %*% par[c(1:2)])
  kappa_param = ifelse(durations < kappa_par[4],
                       exp((cbind(rep(1, n_covariates) , log(durations)) %*% kappa_par[c(1:2)])),
                       exp((cbind(rep(1, n_covariates) , log(durations)) %*% kappa_par[c(5,3)])))
  #sigma
  sigma_param =ifelse(durations < sigma_par[4],
                      exp((cbind(rep(1, n_covariates) , log(durations)) %*% sigma_par[c(1:2)])),
                      exp((cbind(rep(1, n_covariates) , log(durations)) %*% sigma_par[c(5,3)])))

  #xi
  xi_param = ((cbind(rep(1, n_covariates) ,log(durations)) %*% xi_par))
  xi_param[xi_param <= 1e-6] = 1e-6
  #return(list("init" = init, "fitted_params" = par.optim, ,"fits" = station_fit))

  return(list("fitted_params" = par.optim, 'kappa_param' = as.vector(kappa_param), 'scale_param' = as.vector(sigma_param),'shape_param' = as.vector(xi_param)))

}




#' function for fitting of Scaling IDF models
#'
#'  Eight diffent models are considered (simple scaling, General IDF formulation of Koutsoyiannis et al. (1998),
#'  scaling breaks and shape parameter dependence with duration). shape parameter is estimated based on a linear - log  model,
#'   For details, see Haruna et al 2023
#'
#'
#' @param station_data data_frame of observation for a station (nobs x aggregation duration)
#' @param durations  a (vector) of durations (same length as ncol of station_data)
#' @param declustering_duration a vector same length as \code{duration} .Whether the data should be temporally declustered, if yes the time step for each duration, or a vector of 1 otherwise
#' @param initial_params an object returned by \code{egpd_idf_init}. Can be omitted, then the \code{egpd_idf_init} will be implicitly called to find initial values
#' @param censored a scalar or vector of \code{length(durations)}, for the left censoring to be applied to data of each duration. If a scalar is given, it will be divided by \code{durations}
#' @param init_time_step  a scalar, eg 1 or 2. The time step to start declustering the data. eg, for hourly data, if  \code{declustering =3} and  \code{init_time_step = 2}, then the 2nd hour will be selected, and then a sequence is applied
#' @param fitting_method either \code{"mle"} for maximum likelihood or \code{"pwm"} for probability weighted momoments. Defaults to \code{"mle"}
#' @param simple_scaling logical, defaults to \code{TRUE}, defines the initial values values that would be returned
#' @param multi_regime logical, defaults to \code{FALSE}, whether a break in scaling is allowed in the scale parameter
#' @param xi_constant a scalar, defaults to \code{TRUE}, whether the shape parmeter is considered independent of duration or not
#' @param auto_fit logical, whether an automatic fit is to be done to find a censoring threshold that minimizes the normalized root mean square error; defaults to \code{TRUE}
#' @param nrmse_tol a scalar, defaults to \code{0.1}. defaines a relatively good egpd fit. If higher, and \code{auto_fit=TRUE}, automatic fit is done
#' @param simple_scaling logical, defaults to \code{TRUE}, defines the initial values values that would be returned
#' @param use_r_optim logical, defaults  to \code{FALSE}, if yes, \code{optim} will be used to find the best lower censoring values, otherwise, descrete values will be tested
#' @param nrsme_quantile a number in percentage, or any character eg \code{nrsme_quantile = 90}, the normalized root mean square error will only be computed on excesses of the quaneile.  If a charanter, then a weighted normalized root mean square error will be used.
#' @param optim_algo the \code{optim} algtorthm to use. defaults to \code{"Nelder-Mead" }
#'
#' @details
#'   to be added
#'
#' @return # A list:
#'    init: the initial parameters to use for fitting a egpd IDF model
#'    fits: egdp parameters fit for each duration seperately
#'    regres_summary: summary of the lm fit of the sigma vs duration
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
#'  ## get initial values
#' \dontrun{
#'  initial_params = egpd_idf_init(station_data = station_data,
#'                     durations = durations, fitting_method = "mle",
#'                    declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
#'                   auto_fit = T, nrmse_tol = 0.1,use_r_optim = T, nrsme_quantile = 0)
#'  ## fit the data driven IDF
#'  fitted_idf = fit_egpd_idf_scaling_models(station_data = station_data, durations = durations,
#'                  censored = initial_params$fits$lower_threshold,
#'                  declustering_duration = c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72),
#'                  fitting_method = 'mle',  initial_params = initial_params,
#'               simple_scaling = T, multi_regime = T, xi_constant = F,
#'               init_time_step = 1,auto_fit = F, nrmse_tol = 0.1, use_r_optim = F,
#'               nrsme_quantile = 0, optim_algo = "BFGS")
#'
#'  #check 'optim' params, for convergence etc
#'  fitted_idf$fitted_params
#'
#'  #parameters of the IDF
#'  fitted_idf$fitted_params$par
#'
#'  # fitted egpd parameters for the given durations
#'  kappa_fit =  fitted_idf$kappa_param
#'  sigma_fit = fitted_idf$scale_param
#'  xi_fit = fitted_idf$shape_param
#'
#'  #compute nrmse to check quality of fit
#'  nrmse_d = compute_nrsme(station_data, c(1,2,3,6,10,12, 16,18, 24, 48, 72),
#'  kappa_fit, sigma_fit, xi_fit, init_time_step = 1, q = 0)
#'  nrmse_d
#'
#' # Plot the IDF curves
#' plot_egpdidf_curves(station_data = station_data,  kappa_fit = kappa_fit,
#'  sigma_fit = sigma_fit, xi_fit = xi_fit, durations,
#'  declustering_duration=c(1,2,3,6,10,12, 16,18, 24, 48, 72), npy = 92, init_time_step=1 )
#'  }
#' @export
fit_egpd_idf_scaling_models <- function(station_data, durations, declustering_duration, initial_params,
                                        censored,  fitting_method = "mle", simple_scaling = T, multi_regime =F,
                                        xi_constant = T, init_time_step=1, auto_fit = T,nrmse_tol =0.1,use_r_optim=F,
                                        nrsme_quantile =0, optim_algo = "Nelder-Mead" ){



  if (missing(initial_params)) {
    message("initial fittings")
    # fit the data for each duration seperately
    if (use_r_optim) { # r optim function used to find the best lower_censoring; slower but find the best threshold
      station_fit <- local_fit_IDF_h_par(sample = station_data, fitting_method= fitting_method, auto_fit = auto_fit, nrmse_tol = nrmse_tol,
                                         low_th = c(seq(0,1.2,0.1)), durations = durations, declustering_duration = declustering_duration,
                                         init_time_step=init_time_step, q = nrsme_quantile)

    } else{ # here discrete values of lower threshold are tested, its possible to miss the best value
      station_fit <- local_fit_IDF_h(sample = station_data, fitting_method= fitting_method, auto_fit = auto_fit, nrmse_tol = nrmse_tol,
                                     low_th = c(seq(0,1.2,0.1)), durations = durations, declustering_duration = declustering_duration,
                                     init_time_step=init_time_step, q = nrsme_quantile)

    }

    #extract the fitted parpams
    sigma_vec = station_fit$scale_param
    kappa_vec = station_fit$kappa_param
    xi_vec = station_fit$shape_param

    # estimaete the params of the IDF
    # by definiton, only sigma of EGPD shows scaling
    lmsig <- lm(log(sigma_vec)~log(durations))
    sigma0 = exp(coefficients(lmsig)[1]) # intercept of the fit
    eta0 =  -(coefficients(lmsig)[2]) # -ve slope of the fit

    # kappa and xi should be uniform (indepernedent of duration)
    kappa0 = median(kappa_vec)

    #xi is dependent on duration, we fit a spline relationship
    if (!xi_constant){
      xi_fit = fit_linear_model(initial_params$fits$shape_param, durations, model_type = 'linear_log_SR')
    }



    #if xi is taken as INDependent of duration
    xi0 = median(xi_vec)

    #initial theta
    theta0 = 0

    init = c(kappa0, sigma0, xi0, eta0)
    names(init) = c("kappa0", "sigma0", "xi0", "eta")

    #if not simple scaling, then we need to include theta
    if (!simple_scaling) {
      init[5] = theta0
      names(init)[5] = "theta"
    }
  } else{
    init =  initial_params$init
    if (!simple_scaling) {
      init[5] = 0
      names(init)[5] = "theta"
    }
    xi_vec = initial_params$fits$shape_param[colnames(station_data)]
    if(!xi_constant){
      init_xi = lm((xi_vec)~log(durations))$coefficients
      names(init_xi) = c('xi0', 'xi_eta1')
      #here we want the initial intercep > 0 and -ve slope
      if (init_xi['xi0'] < 0)
        init_xi['xi0'] = -init_xi['xi0']
      if (init_xi['xi_eta1'] > 0)
        init_xi['xi_eta1'] = -init_xi['xi_eta1']

      init[c('xi0', 'xi_eta1')] = init_xi
    }

  }
  scaling_break = NULL
  if (multi_regime) {
    #fit segmented, TWO regimes exist
    xxx=fit_segmented(paramet = initial_params$fits$scale_param[colnames(station_data)], covar = durations)
    xxx['eta1'] = max(0.1, xxx['eta1'])
    xxx['eta2'] = max(0.1, xxx['eta2'])
    init = init[!names(init) %in% c('eta', 'theta')]
    init[names(xxx)[2:5]] = xxx[2:5]
    #names(init[4:7]) <- names(xxx[2:5])
    init[2] = xxx[1]
    names(init)[2] = 'sigma01'
    #names(init) = c('kappa0', 'sigma01', 'xi0', 'eta1', 'eta2', 'sigma02', 'k')
    if (!simple_scaling) {
      init["theta"] = 0
      #names(init)[8] = "theta"
    }
    scaling_break = init['k']
    init = init[! names(init) %in% c('sigma02', 'k')]
  }
  if (length(censored)==1)
    censored = censored/durations

  #preparae tha data into a suitable format ----
  n_cens = purrr::map_dbl(seq_along(declustering_duration), function(d){

    data_dec =  decluster_set_NA(x = station_data[[d]], init_time_step = init_time_step, step = declustering_duration[d]) %>% na.omit
    # censoring_value = init_param %>%  filter(duration == durations[d], area == grid_area[a]) %>%  pull(lower_c)
    length(data_dec[data_dec < censored[d]])
  })


  station_data_censored = purrr::map_dfc(seq_along(declustering_duration), function(d){
    data_dec =  decluster_set_NA(x = station_data[[d]], init_time_step = init_time_step, step = declustering_duration[d])

    # seclect obs that are >= censoring value for duration d
    #censoring_value = init_param %>%  filter(duration == durations[d], area == grid_area[a]) %>%  pull(lower_c)
    data_dec[data_dec < censored[d]] =  NA
    data_dec = data.frame( data_dec ) %>%  setNames(durations[d])
    data_dec
  }) %>% #pivot longer
    tidyr::pivot_longer(cols = 1:length(durations),names_to = "duration", values_to = "precip",
                        names_transform = list(duration = as.numeric) ) %>%  arrange(duration)

  station_data_censored = tidyr::drop_na(station_data_censored)

  censored_data = station_data_censored %>%  group_by(duration) %>% group_split()

  # function for the likelihhod (with Xi dependent on D) to be mqximized taking into acoount lower censoring
  LK.EGPD.idf <- function(init, scaling_break, durations, censored,multi_regime, censored_data,  n_cens){
    #### INPUT
    # init (vector(5)): kappa; sigma (scale), xi (shape), eta, theta (IDF of koutsoyiannis, 1998)
    # station_data: data_frame of observation for a station (nobs x aggreation duration)
    # durations: (vector) aggeragtion durations (same lenght as ncol of station_data)
    #censored: (scalar) censoring threshold (lower) to be applied to all durations ()
    # xi_fit : predicted xi from the smoothing spline regression

    #### OUTPUT
    # Likelihood value


    #init_values
    if (multi_regime) {
      init['sigma02'] = init['sigma01']*scaling_break^(init['eta2'] - init['eta1'])
      init['k'] = scaling_break
    }
    nms = names(init)
    if (any(nms =='theta'))
      init = c(init[-which(names(init) =="theta")], init['theta'])
    #init[] = init[c()]



    if (xi_constant) {
      ## ***** check all params are within range
      if (any(nms == 'theta') &  init['theta'] < 0 ) return(1e100) # check that theta remains positive

      if (any(nms == 'eta2') )
        if(init['eta2'] <= 0 | init['eta2']  >= 1) return(1e100) #


      # check that slope eta is within range, named eta or eta1 depending on no of regimes
      if(multi_regime){
        if ( init[c('eta1')] >= 1 | init[c('eta1')] <= 10^(-6)) return(1e100)
      } else {
        if ( init[c('eta')] >= 1 | init[c('eta')] <= 10^(-6)) return(1e100)
      }
      ### ******** end check

      #here only sigma is a function of duration

      # collate all data in single vectors, arranging them by durations eg
      # vec_precip= c(1,1,1..., 2,2...) contains obs form d =1, to d =n
      if (multi_regime) {
        vec_sigma_d= sapply(durations, function(d){
          if (d < init['k']) {
            sigma_d = ifelse(any(nms =='theta'),init['sigma01']/(d+init['theta'])^init['eta1'],init['sigma01']/(d)^init['eta1'])
            #xi_d = medi
          } else{
            sigma_d = init['sigma02']/(d)^init['eta2'] #second regime should remain SS, no theta
          }
          sigma_d
        })


      } else {
        #here  ONE scaling regime exist
        vec_sigma_d= sapply(durations, function(d){
          ifelse(any(nms =='theta'),init['sigma0']/(d+init['theta'])^init['eta'],init['sigma0']/(d)^init['eta'])
        })

      }

      vec_sigma_d_full = sapply(seq_along(durations), function(i) rep(vec_sigma_d[i], length(censored_data[[i]]$precip)) )%>% unlist

      nobs = length(vec_sigma_d_full)
      kappa_vec = rep(init['kappa0'], nobs)
      xi_vec =  rep(init['xi0'], nobs)


      #chek all kappa, xi and all sigma are positive
      if(all(unique(vec_sigma_d)> 0) & init[c('kappa0')] > 0 &
         init[c('xi0')] < 1 & init[c('xi0')] > 10^(-6) ){
        # for left censored data
        contrib.cens1 <- ifelse(any(censored > 0), sum(n_cens[censored>0] * log(p_egpd(x = censored[censored>0], kappa = init['kappa0'], sigma_d  = vec_sigma_d[censored>0], xi = init['xi0']))), 0)
        #for uncensored data
        #contrib.not.cens <- sum(d_egpd(x = vec_p, kappa = kappa_vec, sigma_d = vec_sigma_d, xi = xi_vec))
        contrib.not.cens <- sum(d_egpd(x = sapply(seq_along(durations), function(i) censored_data[[i]]$precip) %>% unlist,
                                       kappa = kappa_vec, sigma_d = vec_sigma_d_full, xi = xi_vec))
        #sum the two likelihhods

        return(-(contrib.cens1+ contrib.not.cens))
      }else{
        return(1e100)
      }
    } else{ ###### XI a function ot duration

      ##  ****** checck that params are within range
      #if (init['xi0'] < 0 | init['xi_eta1'] > 0) return(1e100)
      # if (length(init) == 4 &  init[4] < 0) return(1e100) # check that theta remains positive
      if (any(nms == 'theta') &  init['theta'] < 0 ) return(1e100) # check that theta remains positive

      if (any(nms == 'eta2') )
        if(init['eta2'] <= 0 | init['eta2']  >= 1) return(1e100)

      # check that slope eta is within range, named eta or eta1 depending on no of regimes
      if(multi_regime){
        if (any(init[c('kappa0', 'eta1')] <= 0) | init[c('eta1')] >= 1 | init[c('eta1')] <= 10^(-6) ) return(1e100)
      } else {
        if (init[c('eta')] >= 1 |  init[c('eta')] <= 10^(-6) |init[c('kappa0')] <= 0  ) return(1e100)
      }
      # ******** end check


      # collate all data in single vectors, arranging them by durations eg
      # vec_precip= c(1,1,1..., 2,2...) contains obs form d =1, to d =n
      if (multi_regime) {
        vec_sigma_d= sapply(durations, function(d){
          if (d < init['k']) {
            sigma_d = ifelse(any(nms =='theta'),init['sigma01']/(d+init['theta'])^init['eta1'],init['sigma01']/(d)^init['eta1'])
            #xi_d = medi
          } else{
            sigma_d = init['sigma02']/(d)^init['eta2'] #second regime should remain SS, no theta
          }
          sigma_d
        })
      } else {
        #here  ONE scaling regime exist
        vec_sigma_d= sapply(durations, function(d){
          ifelse(any(nms =='theta'),init['sigma0']/(d+init['theta'])^init['eta'],init['sigma0']/(d)^init['eta'])
        })


      }
      vec_xi_d= sapply(durations, function(d){c(1 ,log(d)) %*% init[c('xi0','xi_eta1')]})
      vec_xi_d[vec_xi_d<1e-6] = 1e-6


      vec_sigma_d_full = sapply(seq_along(durations), function(i) rep(vec_sigma_d[i], length(censored_data[[i]]$precip)) )%>% unlist
      vec_xi_d_full = sapply(seq_along(durations), function(i) rep(vec_xi_d[i], length(censored_data[[i]]$precip)) )%>% unlist

      nobs = length(vec_sigma_d_full)
      kappa_vec = rep(init['kappa0'], nobs)


      # do a likelihood fit (with left censoring, depending on whether the censored value is set to zero or otherwise)




      if(all(unique(vec_sigma_d)> 0)  & all(unique(vec_xi_d)>= 1e-6) & all(unique(vec_xi_d) < 0.5) ){
        # for left censored data
        contrib.cens1 <- ifelse(any(censored > 0), sum(n_cens[censored>0] * log(p_egpd(x = censored[censored>0], kappa = init['kappa0'], sigma_d  = vec_sigma_d[censored>0], xi = vec_xi_d[censored>0]))), 0)
        #for uncensored data
        contrib.not.cens <- sum(d_egpd(x = sapply(seq_along(durations), function(i) censored_data[[i]]$precip) %>% unlist,
                                       kappa = kappa_vec, sigma_d = vec_sigma_d_full, xi = vec_xi_d_full))
        #sum the two likelihhods
        return(-(contrib.cens1+ contrib.not.cens))
      }else{
        return(1e100)
      }
    }

  }
  # optimize the function above
  #message("IDF fitting")
  par.optim = optim(par=init,fn=LK.EGPD.idf,gr=NULL, scaling_break = scaling_break, multi_regime = multi_regime, durations = durations, censored= censored,
                    censored_data =censored_data,  n_cens = n_cens,  control = list(maxit = 3000), hessian = FALSE,method=optim_algo)

  #update the optime parameters to retuen the scaling break and sigma02
  if (multi_regime) {
    par.optim$par['sigma02'] = par.optim$par['sigma01']*scaling_break^(par.optim$par['eta2'] - par.optim$par['eta1'])
    par.optim$par['k'] = scaling_break
    if (any(names(par.optim$par) =='theta'))
      par.optim$par = c(par.optim$par[-which(names(par.optim$par) =="theta")], par.optim$par['theta'])
  }

  #esimate parameters
  par = par.optim$par
  sigma_fit =c()
  kappa_fit = rep(par["kappa0"], length(durations))

  if(xi_constant) {
    xi_fit = rep(par["xi0"], length(durations))
  }else {
    n_covariates = length(durations)
    #xi
    xi_fit = ((cbind(rep(1, n_covariates) ,log(durations)) %*% par[c('xi0','xi_eta1')]))
    xi_fit[xi_fit <= 1e-6] = 1e-6
  }

  for (id in  1:length(durations)) {
    d =  durations[id]
    sigma_fit[id] = get_sigma_d(par = par, d = d)
  }

  return(list("init" = init, "fitted_params" = par.optim, xi_spline = xi_fit, "scale_param" = sigma_fit, "shape_param" = as.vector(xi_fit), "kappa_param" = kappa_fit))

}


fit_segmented = function(paramet, covar){
  yy = log(paramet)
  xx =log(covar)
  lm_fit = lm(yy~xx)


  zz = segmented(lm_fit,   control = seg.control(display = FALSE, alpha = 0.1), psi = log(6))
  #breakpoint
  k = exp(zz$psi[1,2])
  #slopes of line 1 and 2
  eta1 = -slope(zz)$x[1,1]
  eta2 = -slope(zz)$x[2,1]

  #intercept of line 1 and 2
  #sigma02 = sigma01*k^-(eta1-eta2)
  sigma01 = exp(intercept(zz)$x[1,1])
  sigma02 = exp(intercept(zz)$x[2,1])
  result = c(sigma01, eta1, eta2, sigma02, k)
  names(result) = c('sigma01', 'eta1', 'eta2', 'sigma02', 'k')
  return(result)
}
#xxx=fit_segemeneted(paramet = sigma_mat["ZER",-c(1:2),1], covar = durations)



#comparison
#q is specified in two ways:
#in number % eg 0% for all postive data, and 95% for exceedances of a 95% quantile. Quantile is over ALL  values (zeros included)
# as a charcter, for a weighted NMSE, here the wights are inverse of the CDF. (Beguería and Vicente-Serrano [2006]: panthou et al 2012)
nrmse_idf= function(data, kappa, sigma, xi, ranking_method = "average", q = 0){
  data = na.omit(data)
  data = sort(data[data>0], decreasing = F)
  rank_i = rank(data, ties.method = ranking_method)
  q_emp  <- data.frame(d =  data, probs = rank_i/(length(data[data>0])+1))
  q_theo <- qextgp(p=q_emp$probs, prob=NA, kappa=kappa, delta=NA, sigma= sigma, xi=xi , type=1)
  if (is.numeric(q)) {
    threshold = quantile(data,q/100, na.rm = T) #threshold over all values including zeros
    ex_id = which(q_emp$d > threshold) # fidn the index of the excecedences of q% quantile
    rmse <-  sqrt(mean((q_theo[ex_id] - q_emp$d[ex_id])^2))/mean(q_emp$d[ex_id], na.rm=T)
    return(rmse)
  }else{
    weight =  1/(1-q_emp$probs)
    rmse <-  sqrt(mean((q_theo*weight - q_emp$d*weight)^2))/mean(q_emp$d*weight, na.rm=T)
    return(rmse)
  }
}
#function for LINEAR modlels

get_vector_of_param = function(n ,d, covariate_d, init_param){
  #sigma
  if (d < init_param[4]) {
    param_d = ((cbind(rep(1, n) ,log(covariate_d)) %*% init_param[c(1,2)]))
  } else {
    param_d = ((cbind(rep(1, n) ,log(covariate_d)) %*% init_param[c(5,3)]))
  }
  param_d
}

get_vector_of_param2 = function( d, init_param, param  = ""){
  #sigma
  if (d < init_param[4]) {
    param_d = c(1,log(d)) %*% init_param[c(1,2)]
  } else {
    param_d =  c(1,log(d)) %*% init_param[c(5,3)]
  }
  if(param == "xi"){
    return(as.vector(param_d))
  } else{
    return(as.vector(exp(param_d)))
  }

}
get_second_intercept = function(init_param){
  #init_p['sigma01'] + log(init_p['k'])*((init_p['eta1']-init_p['eta2']))
  # init_param[1] + log(init_param[4])*(init_param[3]-init_param[2])
  param02 = init_param[1] + log(init_param[4])*(init_param[2]-init_param[3])
  #names(param02) = sub('01','02',names(param02) )
  param02
  #init_param[1]*(init_param[4])^(init_param[3]-init_param[2])
}

find_init_linear_models = function(paramet, covar, model_type= 'log_log_MR', param_name){
  init_breakponts = log(6)
  if (model_type == 'log_log_MR'){
    m4 = lm(log(paramet)~log(covar))
    y = log(paramet)
    x =log(covar)
    fit = m4
    try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts))
    k = exp(fit$psi[1,2])
    #slopes of line 1 and 2
    eta1 =   slope(fit)$x[1,1]
    eta2 =   slope(fit)$x[2,1]

    if (param_name == "sigma") {
      eta1 =  min(-0.1, slope(fit)$x[1,1])
      eta2 =  min(-0.1, slope(fit)$x[2,1])
    }
    #intercept of line 1 and 2
    #sigma02 = sigma01*k^-(eta1-eta2)
    sigma01 = (intercept(fit)$x[1,1])
    #sigma02 = exp(intercept(fit)$x[2,1])
    result = c(sigma01, eta1, eta2,  k)
    names(result) = paste0(param_name,  c('01', '_eta1', '_eta2',  '_k'))
  } else if (model_type == 'linear_log_MR'){
    m4 = lm(paramet~log(covar))
    y = paramet
    x =log(covar)
    fit = m4
    try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts))
    k = exp(fit$psi[1,2])
    #slopes of line 1 and 2
    eta1 =   slope(fit)$x[1,1]
    eta2 =   slope(fit)$x[2,1]

    if (param_name == "sigma") {
      eta1 =  min(-0.1, slope(fit)$x[1,1])
      eta2 =  min(-0.1, slope(fit)$x[2,1])
    }
    #intercept of line 1 and 2
    #sigma02 = sigma01*k^-(eta1-eta2)
    sigma01 = (intercept(fit)$x[1,1])
    #sigma02 = exp(intercept(fit)$x[2,1])
    result = c(sigma01, eta1, eta2,  k)
    names(result) = paste0(param_name,  c('01', '_eta1', '_eta2',  '_k'))
  }

  return(result)
}


fit_linear_model = function(paramet, covar, model_type){
  init_breakponts = log(6)

  if (model_type == 'linear_log_SR') {
    fit = lm((paramet)~log(covar))
    predicted = fitted.values(fit)
  } else if (model_type == 'log_log_SR'){
    fit = lm(log(paramet)~log(covar))
    predicted = exp(fitted.values(fit))
  }else if (model_type == 'linear_SR'){
    fit = lm((paramet)~(covar))
    predicted = (fitted.values(fit))
  }else if (model_type == 'linear_MR'){
    m4 = lm((paramet)~(covar))
    y = (paramet)
    x = (covar)
    fit = m4
    try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = exp(init_breakponts)), silent = T)
    predicted = (fitted.values(fit))
  }else if (model_type == 'log_log_MR'){
    m4 = lm(log(paramet)~log(covar))
    y = log(paramet)
    x =log(covar)
    fit = m4
    try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    predicted = exp(fitted.values(fit))
  }  else if (model_type == 'linear_log_MR'){
    m4 = lm((paramet)~log(covar))
    y = (paramet)
    x =log(covar)
    fit = m4
    try(fit <- segmented(lm(y~x),   control = seg.control(display = FALSE, alpha = 0.1), psi = (init_breakponts)), silent = T)
    predicted = (fitted.values(fit))
  }
  return(predicted)

}


#function for rolling mean
roll=function (data, ds, which.stations = NULL, names = c("date", "RR"))
{

  agg.station <- function(station) {
    data.s <- data[[station]]
    if (!is.data.frame(data.s)) {
      stop("Elements of 'data' must be data.frames. But element ",
           station, " contains: ", class(data.s))
    }
    names = colnames(data.s)
    if (sum(is.element(names[1:2], names(data.s))) != 2) {
      stop("Dataframe of station ", station, " does not contain $",
           names[1], " or $", names[2], ".")
    }
    dtime <- as.numeric((data.s[, names[1]][2] - data.s[,  names[1]][1]))
    agg.ts <- function(ds) {
      runsum = RcppRoll::roll_sum(data.s[, names[2]], round(ds/dtime),
                                  fill = NA, align = "right")
      runsum <- runsum/(ds/dtime)
      return(runsum)
    }
    data.agg <- lapply(ds, agg.ts)
    df <- do.call(rbind, data.agg)
    return(df)
  }
  if (is.null(which.stations)) {
    which.stations <- if (is.null(names(data))) {
      1:length(data)
    }
    else {
      names(data)
    }
  }

  station.list <- lapply(which.stations, agg.station)
  names(station.list) = names(data)
  stationdf = t(do.call("rbind", station.list))
  return(data.frame(stationdf))
}
