#'  Function to fit egpdIDAF model
#'
#'  A regression based (data driven) IDAF model based on the egpd distribution.
#'  IDAF stands for Intensity-Duration-Area-Frequency.
#'
#' @param st_data_ad a list of size \code{length(grid_area)}. Each containing 1) a df of the
#' intensities (of various durations) of the given aggregation area, and
#' b) a vector of size \code{length(durations)}  containing the number of observations below a given censoring threshold.
#' @param durations  a (vector) of durations
#' @param grid_area  a (vector) of aggregation areas
#' @param formula  a list (size = 3) of  formulars, one each for each egpd parameter
#' @param use_mle_init logical, defaults  to \code{FALSE}, if yes, an iterative pairwise likelihood fitting is done. See ...
#'
#' @details
#'   to be added
#'
#' @return # A list:
#'
#'    fits: egdp parameters fit for each duration seperately
#'    optim details, to check convergence
#'
#' @examples
#'  ## load the data
#'  data(precip_areal)
#'   side_vec = c(0,1,2,3,4,6,8,11,13,16)
#'   grid_area = (1+2*side_vec)^2
#'   area_names = paste0("A_",grid_area,"km2")

#'  durations = c( 1,2, 3,  6,  10, 12,  16,  24, 48, 72)
#'  duration_names= c(paste0("D_", durations, "h"))
#'  declustering_duration = c( 3, 4, 5,  8,  10, 12,  16,   24, 48, 72)
#'
#'  #extract the data  at the pixel location
#'
#' \dontrun{
#'
#'  #1.0 obtaining initial values. This is done by fitting EGPD locally to each data of A and D. ------
#'  ## Read  the 'help' of 'egpd_idf_init' function in the  'egpdIDF' package .
#'
#'  # craere an emty array to store the parameters
#'   param_matrix = array(data = NA, dim = c( 5, length(durations), length(side_vec)),
#'                     dimnames = list( c("kappa", "sigma", "xi", "nrmse", "lower_c"),
#'                                      durations, grid_area))
#'  #iterate over the areal data, each time  aggregate the data into various durations, and fit egpd
#'  for (a in seq_along(area_names)) {
#'    cat(-a)
#'    station_data= aggregate_data(sample_data = areal_data,
#'          st_code = colnames(areal_data)[a+1], durations)
#'
#'    #--- fit the egod:  Read  the 'help' of 'egpd_idf_init' function in the  'egpdIDF' package
#'    initial_params = egpd_idf_init(station_data = station_data,
#'                                          durations = durations,
#'                                          fitting_method = "mle",
#'                                          declustering_duration =  declustering_duration,
#'                                          auto_fit = F, use_r_optim = F)
#'    #  store the params
#'    param_matrix[, ,a] = rbind(initial_params$fits$kappa_param, initial_params$fits$scale_param,
#'                             initial_params$fits$shape_param, initial_params$fits$nrsme,
#'                             initial_params$fits$lower_threshold)
#'  }
#'
#'  # 2. fitting the EGPD-IDAF  model -------
#'  ## 2.1  putting the data in the right format ----
#'  init_param = as.data.frame.table(param_matrix[1:5,,]) %>%
#'    set_names(c("param", "duration", "area", "value") )  %>%
#'     pivot_wider(names_from = "param",values_from = "value") %>%
#'    mutate_if(is.factor, ~ as.numeric(as.character(.)))  %>%  arrange(duration)
#'
#'  st_data_ad =  map(seq_along(area_names), function(a) {
#'    # Aggregate the data to durations
#'    xx <-  aggregate_data(sample_data = areal_data,
#'       st_code = colnames(areal_data)[a+1], durations)
#'    names(xx) <- durations
#'
#'    # Decluster each data and count the number of censored observations
#'    n_cens <- sapply(seq_along(declustering_duration), function(d) {
#'      data_dec =  decluster_set_NA(x = xx[[d]], init_time_step = 1,
#'       step = declustering_duration[d]) %>% na.omit
#'     #data_dec <- decluster_set_NA(xx[[d]], init_time_step, step = d)
#'     censoring_value <- init_param$lower_c[init_param$duration == durations[d] &
#'       init_param$area == grid_area[a]]
#'      sum(data_dec > 0 & data_dec < censoring_value)
#'    })
#'
#'    # Filter and pivot the declustered data
#'    xx1 <- map(seq_along(declustering_duration), function(d) {
#'      data_dec <- decluster_set_NA(xx[[d]], init_time_step =1, step = declustering_duration[d])
#'      censoring_value <- init_param$lower_c[init_param$duration == durations[d]
#'       & init_param$area == grid_area[a]]
#'     data_dec[data_dec < censoring_value] <- NA
#'      data.frame(precip = data_dec, duration = durations[d], area = grid_area[a])
#'    }) %>%  do.call(what = rbind) %>%  arrange(duration) %>% drop_na
#'
#'    # Return the filtered and pivoted data and the number of censored observations
#'    list(data = xx1, ncens = n_cens)
#'})
#'
#'  ## 2.2  specifying the relationship ----------
#'  kappa_formula = "log(kappa) ~ area + log(area)*sqrt(duration)  + log(area)*duration"
#'  sigma_formula = "log(sigma) ~ area + duration*log(area) + log(area)*sqrt(duration)"
#'  xi_formula =  " xi ~ log(area)*sqrt(duration) + log(area)*duration"
#'
#'  formula = list(kappa_formula, sigma_formula,  xi_formula)
#'
#'  ## 2.3 fitting the model ---------
#'  fitted_idaf = fit_idaf_dd(st_data_ad = st_data_ad, init_param = init_param,
#'       durations = durations, grid_area = grid_area,formula = formula, use_mle_init = F)
#'
#'  kappa_mat_fitted = fitted_idaf$kappa_param
#'  sigma_mat_fitted = fitted_idaf$scale_param
#'  xi_mat_fitted = fitted_idaf$shape_param
#'
#'  ## 2.4  diagonostics (qqplots) of the model ----------
#'  for (a in seq_along(area_names)) {
#'    cat(-a)
#'   station_data= aggregate_data(sample_data = areal_data,
#'       st_code = colnames(areal_data)[a+1], durations)
#'   pp = map(seq_along(durations), function(d){
#'          qqplot_egpd(durations[d]*station_data[seq(1, nrow(station_data),
#'                          declustering_duration[d]),d],
#'                          kappa= kappa_mat_fitted[[a+1]][d],
#'                          sigma = sigma_mat_fitted[[a+1]][d]*durations[d],
#'                           xi = xi_mat_fitted[[a+1]][d])
#'      })
#'    p= ggarrange(plotlist = pp, ncol = 3, nrow = 4)
#'    print(p)
#'
#'  }
#' #}

#' @export

fit_idaf_dd =  function(st_data_ad, init_param,  durations, grid_area,formula, use_mle_init = F){
  #st_data_ad, : a list containing the data in a particular format
  #init_param, ; a df of initial paramters
  #durations, : a vector of aggregation durations
  #grid_area,: a vector of area
  #formula,  : a list of s3 formulars, one each for an egpd paramter
  #use_mle_init: whether  a pairwise likelihhod should be done to estimate the params

  if (missing(init_param)) {
    stop("arg init_param missing")
  }
  # prelimnaries  ------
  #censored data
  censored_data = map_dfr(seq_along(grid_area), function(i)  st_data_ad[[i]]$data ) %>%
    group_by(duration, area) %>% group_split()
  #number of censosred data
  n_censored =  map_dfc( seq_along(grid_area), function(i) {
    x=  st_data_ad[[i]]$ncens
    data.frame(x) %>% setNames(grid_area[i])
  })  %>% mutate(duration = durations) %>%
    tidyr::pivot_longer(cols = 1:length(grid_area),names_to = "area", values_to = "n_cens",
                        names_transform = list(area = as.numeric) )

  # lm fits
  lm_fits = map(1:3, function(p){
    lm(formula[[p]], data = init_param, x = T)
  })

  # covariate matrix (X matrix)
  cov_mat = map(1:3, function(p){
    x_mat = lm_fits[[p]]$x
  }) %>% set_names(c("kappa","sigma", "xi"))

  #init values of the regression coefficients
  init = map(1:3, function(p){
    pname = c("kappa","sigma", "xi")[p]
    init_val = lm_fits[[p]] %>% coef
    names(init_val) = paste0(substr(pname,1,2),"_", names(init_val))
    init_val
  }) %>% flatten_dbl()

  # link
  link_type = map_chr(1:3, function(p){
    link = "identity"
    is.log = any(lm_fits[[p]]$terms[[2]] %>% as.character() == "log")
    if (is.log) {
      link = "log"
    }
    link
  }) %>% set_names(c("kappa","sigma", "xi"))

  # likelhood function  ------
  LK_EGPD_idaf_dd <- function(init, free_params, cov_mat, censored_data, n_censored, init_param, link_type){

    free_params[names(init)] = init

    kappa_init = free_params[which(substr(names(free_params), 1, 2) == "ka")]
    sigma_init = free_params[which(substr(names(free_params), 1, 2) == "si")]
    xi_init = free_params[which(substr(names(free_params), 1, 2) == "xi")]

    vec_kappa_d = get_fitted_param(x_mat = cov_mat[["kappa"]], beta_vec = kappa_init, link = link_type["kappa"])
    vec_sigma_d = get_fitted_param(x_mat = cov_mat[["sigma"]], beta_vec = sigma_init, link = link_type["sigma"])
    vec_xi_d = get_fitted_param(x_mat = cov_mat[["xi"]], beta_vec = xi_init, link = link_type["xi"])

    if (all(vec_xi_d < 0) |  any(vec_xi_d > 0.5) ) {
      return(1e100)
    }else{
      vec_xi_d[vec_xi_d<1e-6] = 1e-6
    }
    #vec_xi_d[vec_xi_d<1e-6] = 1e-6
    #vec_kappa_d[vec_kappa_d<0.1] = 0.1
    #vec_sigma_d[vec_sigma_d<0.1] = 0.1

    #vec_xi_d[vec_xi_d<1e-6] = 1e-6
    # likelhood
    if(all(c(vec_kappa_d, vec_sigma_d)> 0)){
      is_censored = init_param$lower_c >0
      contrib.cens1 <- ifelse(any(init_param$lower_c > 0), sum(n_censored$n_cens[is_censored] * log(p_egpd(x = init_param$lower_c[is_censored], kappa = vec_kappa_d[is_censored],
                                                                                                           sigma_d  = vec_sigma_d[is_censored], xi = vec_xi_d[is_censored]))), 0)

      #for uncensored data
      contrib.not.cens <- map_dbl(seq_along(vec_sigma_d), function(i){
        sum(d_egpd(x = censored_data[[i]]$precip, kappa = vec_kappa_d[i], sigma_d = vec_sigma_d[i], xi = vec_xi_d[i]), na.rm = T)
      }) %>%  sum

      #print(c(contrib.cens1, contrib.not.cens))
      #print(c(min(vec_kappa_d), min(vec_sigma_d), min(vec_xi_d)))
      #print(init)
      #sum the two likelihhods
      nlklhood_sum = -(contrib.cens1+ contrib.not.cens)
      #check if its finite
      if(is.nan(nlklhood_sum) | is.infinite(nlklhood_sum)) {
        return(1e100)
      } else{
        return(nlklhood_sum)
      }

    } else{
      return(1e100)
    }


  }

  #optimization ----------
  if (use_mle_init) {

    # par.optim =  optim(par=init,fn=LK_EGPD_idaf_dd,gr=NULL, free_params = init, cov_mat= cov_mat, censored_data = censored_data, n_censored = n_censored,
    #                    init_param = init_param, link_type = link_type,
    #                    control = list(maxit = 10000), hessian = FALSE,method="BFGS")
    # init = init
    i = 1
    tol = 100000
    likl = LK_EGPD_idaf_dd(init = init  , cov_mat= cov_mat, free_params = init,censored_data = censored_data, n_censored = n_censored, init_param = init_param, link_type = link_type)

    init_temp = init
    while (tol > 0.05) {
      #message(i)
      for (pp in c("ka", "si", "xi")) {
        pname =  which(substr(names(init), 1, 2) == pp)
        par.optim = optim(par=init_temp[pname],fn=LK_EGPD_idaf_dd,gr=NULL, free_params = init_temp, cov_mat= cov_mat, censored_data = censored_data, n_censored = n_censored,
                          init_param = init_param, link_type = link_type,
                          control = list(maxit = 10000), hessian = FALSE,method="BFGS")
        init_temp[pname] = par.optim$par
        #print(c(par.optim$value, par.optim$convergence))
      }
      tol = likl -  par.optim$value
      print(tol)
      likl = par.optim$value
      i =i +1
    }
    init_ml = init_temp
    par.optim =  optim(par=init_ml,fn=LK_EGPD_idaf_dd,gr=NULL, free_params = init_ml, cov_mat= cov_mat, censored_data = censored_data, n_censored = n_censored,
                       init_param = init_param, link_type = link_type,
                       control = list(maxit = 10000), hessian = FALSE,method="BFGS")


  } else {
    par.optim =  optim(par=init,fn=LK_EGPD_idaf_dd,gr=NULL, free_params = init, cov_mat= cov_mat, censored_data = censored_data, n_censored = n_censored,
                       init_param = init_param, link_type = link_type,#lower = init*0.1, upper = init*1.9,
                       control = list(maxit = 10000), hessian = FALSE,method="BFGS")
  }
  #LK_EGPD_idaf_dd(init = init_ml  , cov_mat= cov_mat, free_params = init_ml,censored_data = censored_data, n_censored = n_censored, init_param = init_param, link_type = link_type)

  # prepare output ------

  par = par.optim$par
  #prdictions
  kappa_par = par[which(substr(names(par), 1, 2) == "ka")]
  sigma_par = par[which(substr(names(par), 1, 2) == "si")]
  xi_par = par[which(substr(names(par), 1, 2) == "xi")]


  kappa_param = get_fitted_param(x_mat = cov_mat[["kappa"]], beta_vec = kappa_par, link = link_type["kappa"])  %>% get_param_df_wide
  sigma_param = get_fitted_param(x_mat = cov_mat[["sigma"]], beta_vec = sigma_par, link = link_type["sigma"])   %>% get_param_df_wide
  xi_param = get_fitted_param(x_mat = cov_mat[["xi"]], beta_vec = xi_par, link = link_type["xi"])  %>% get_param_df_wide
  xi_param[xi_param<1e-6] = 1e-6
  #kappa_param = cov_mat[["kappa"]] %*% kappa_par %>% as.vector %>% get_param_df_wide
  #sigma_param = cov_mat[["sigma"]]%*% sigma_par %>% as.vector  %>% get_param_df_wide
  #xi_param =  cov_mat[["xi"]] %*% xi_par %>% as.vector  %>% get_param_df_wide




  return(list("fitted_params" = par.optim, 'kappa_param' = kappa_param, 'scale_param' = sigma_param,'shape_param' = xi_param))

}
# init = c(0.713658872, 0.011082547,-0.088001442,-0.148979924,-0.224580097,0.007504015,-0.312668515,0.057917281,0.073191910,-0.144571914,0.154911257,-0.014873820,  0.015953700 )
# LK_EGPD_idaf_dd(init = init  , cov_mat= cov_mat, free_params = init,censored_data = censored_data, n_censored = n_censored, init_param = init_param, link_type = link_type)
# par.optim =  nloptr(x0 = init,eval_f = LK_EGPD_idaf_dd, free_params = init, cov_mat= cov_mat, censored_data = censored_data, n_censored = n_censored,
#                    init_param = init_param, link_type = link_type, opts = list("algorithm"="NLOPT_LN_COBYLA", "maxeval" =1000))


#' @export
#functions ----
#decluster a time series temporally, by replacing all the skipped time steps with NA
decluster_set_NA <- function(x, init_time_step, step) {
  # Get indices of valid and invalid time steps
  val_id <- seq(init_time_step, length(x), step)
  inval_id <- setdiff(seq_along(x), val_id)

  # Replace invalid time steps with NA
  x[inval_id] <- NA

  # Replace non-positive values with NA
  x[x <= 0] <- NA

  return(x)
}


#' @export
#convert a datafrane to long form
get_param_df_long  = function(param_matrix,  new_name = "", grid_area_name =grid_area , durations_names = durations){
  if(new_name == "") new_name = "param"
  df = param_matrix %>% data.frame() %>% setNames(grid_area_name) %>%
    mutate( duration = durations_names) %>%
    tidyr:: pivot_longer(cols = seq_along(grid_area_name), names_to ="area", values_to = new_name, names_transform = list(area = as.numeric))
  return(df)
}


#' @export
#convert a datafrane to wide form
get_param_df_wide =  function(param_val){
  df = tibble(duration = rep(durations,  each =length(grid_area)), area =  rep(grid_area,  length(durations)), "param" = param_val) %>%
    tidyr::pivot_wider(names_from = "area", values_from = "param" )
  df

}


#' @export
get_fitted_param = function(x_mat, beta_vec, fitted = NA, link = "identity"){
  if (missing(x_mat) | missing(beta_vec)) {
    y_fitted =  fitted
  } else{
    y_fitted = x_mat %*% beta_vec %>% as.vector
  }

  if (link == "log")
    y_fitted = exp(y_fitted)

  return(y_fitted)
}


#' @export
#apliyng the selection
get_link_type = function(model){
  link = "identity"
  is.log = any(model$terms[[2]] %>% as.character() == "log")
  if (is.log) {
    link = "log"
  }
  link
}
