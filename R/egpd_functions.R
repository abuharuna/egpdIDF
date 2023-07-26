#' @export
#returns  the negative loglikelohood (log of the density function) of egpd for sigma = sigma_D
d_egpd=function(x, sigma_d, xi, kappa){
  if (all(unique(xi) > 0.000001)) {
    #for xi greater than zero
    -(log(sigma_d/kappa)+ (1+1/xi)*log(1+xi*x/sigma_d) - (kappa-1)*log(1-(1+xi*x/sigma_d)^(-1/xi)) )
  }else{
    # for xi = 0
    -(log(sigma_d/kappa) + x/sigma_d - (kappa-1)*log(1-exp(-x/sigma_d)))
  }
}
# returens the PDF of egpd
p_egpd=function(x, sigma_d, xi, kappa){
  if (all(unique(xi) > 0.000001)) {
    #for xi greater than zero
    (1-(1+xi*x/sigma_d)^(-1/xi))^(kappa)
  }else{
    # for xi = 0
    (1-exp(-x/sigma_d))^kappa
  }
}

#' @export
# Probability density function of duration-dependent EGPD distribution
degpd_d <- function(x,  kappa, sigma0, xi, eta, theta, d) {

  if (any(c(length(kappa), length(sigma0), length(xi), length(theta),
            length(eta), length(d)) > 1)) {
    message("One of the parameters mut, sigma0, xi, theta, eta, tau is a vector. ",
            "This is not intended and might cause an error.")
  }
  if (any(d + theta <= 0)) {
    warning("Some shape parameters are negative,resulting from a negativ theta ",
            theta, " this will prododuce NAs.")
  }
  if (d + theta <= 0) {
    return(rep(NA, length(q)))
  }
  else {
    sigma_d =  sigma0/(d + theta)^eta
    dd = dextgp(x = x, kappa = kappa, sigma = sigma_d, xi = xi, type = 1, log = F)
    return(dd)

  }

}


#function to extract the fitted sigma_d from a scaling model
get_sigma_d <- function(par, d){
  #find if single or multiple regime based on the number of parametres
  nms = names(par)
  multi_regime = ifelse(any(nms == 'eta1'),   T , F )
  #find if xi varies with duration based on the existance of the "xi0"
  # xi_constant = ifelse(any(nms == "xi0"), T, F)

  if (multi_regime) {

    if (d <= par['k']) {
      # check to see if theta esixt or not, if yes then koutsiyyanis model
      sigma_d = ifelse(any(nms == 'theta'),par['sigma01']/(d+par['theta'])^par['eta1'], par['sigma01']/(d)^par['eta1'])
    } else{
      sigma_d = par['sigma02']/(d)^par['eta2']
    }

  } else{ #single regime
    sigma_d = ifelse(any(nms == 'theta'),par['sigma0']/(d+par['theta'])^par['eta'], par['sigma0']/(d)^par['eta'])
    # if (xi_constant) {
    #   sigma_d =  ifelse(length(par)==4,par[2]/(d)^par[4] ,par[2]/(d+par[5])^par[4])
    # } else { # xi is not constant
    #   sigma_d =  ifelse(length(par)==3,par[2]/(d)^par[3] ,par[2]/(d+par[4])^par[3])
    # }
  }
  return(sigma_d)
}

#' @export
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
