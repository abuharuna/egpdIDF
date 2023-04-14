# egpdIDF: Modeling Intensity-Duration-Frequency (IDF) Curves using the Extended Genaralized Pareto Distribution (EGPD)

# Introduction

Intensity-Duration-Frequency (IDF) curves provide the link between precipitation  intensity, duration, and non-exceedance frequency (or rather the return period). It is a  very common and useful tool in the area of water resources engineering. IDF curves are  practically used to infer high return levels of rainfall intensities for the hydrological designs of structures such as sewer lines, culverts, drains, dams, dykes, etc. 

IDF curves are usually modeled using Generalized Extreme Value distributions (GEV). Here, we consider the three parameter Extended Generalized Pareto Distribution (EGPD) of Naveau et al 2016 to build IDF curves.  EGPD models the distribution of the non-zero rainfall intensities, not only extremes, thereby  making efficient use of the available information. The package provides the implementation of 10 different IDF modelling approaches as contained in Haruna et al 2023.


# Runnig the Functions

## Loading the package

The package is loaded by calling the function below:

```{r setup}
library(egpdIDF)
```



## Aggeregation precipitation data to intensities of longer durations

Before the IDF curves are modeled, the data has to be aggregated to intensities of different duration that are specified by the user. 

```{r, eval = FALSE}

 ## load the data
 data("precipdata")

 ## Here the resolution of the data is in 'hours', we want to aggeregate the data up to 72 hours
 ## specify the aggregation durations

 durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)

 ## get the aggrageted data for each of the
 station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",  durations = durations)

 head(station_data)
 

}
```


## Fitting EGPD to intensities of each duration separately

We can fit EGPD to the intensities of each duration separately. This helps to investigate how each EGPD parameter varies with duration. We can also investigate the quality of the EGPD fit. Finally, the outputs of the function are used to initialize the other functions ues in modeling the IDF curves.

```{r, eval=FALSE}
initial_params = egpd_idf_init(station_data = station_data,
                durations = durations, fitting_method = "mle",
                declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72), auto_fit = FALSE)
 ## check the fitted egpd parameters for each duration
 initial_params$fits$kappa_param
 initial_params$fits$scale_param
 initial_params$fits$shape_param

 ## check the quality of the fit
 initial_params$fits$nrsme

 ## for a good fit, 'nrsme" should be small.

 ##  Tt's always good to  use left censoring with 'mle' fit. Lets try, and check the 'nrmse' again
 ##  we set "auto_fit=T", "nrmse_tol=0.1"  "use_r_optim=T" , "nrsme_quantile = 0".
 ## Check the arguments for their meaning
 initial_params = egpd_idf_init(station_data = station_data,
           durations = durations, fitting_method = "mle",
            declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
              auto_fit = TRUE, nrmse_tol = 0.1,use_r_optim = TRUE, nrsme_quantile = 0)
 ## check the quality of the fit
 initial_params$fits$nrsme
 ## check the parameters
 initial_params$fits$kappa_param
 initial_params$fits$scale_param
 initial_params$fits$shape_param
```


## Fitting data-driven IDF models

This class of IDF model rely on empirically determined functional relationships between the EGPD parameter and duration. It is not based on any physical model assumption but rather based on the identified empirical relationship.

```{r, eval=FALSE}
## load the data
 data("precipdata")

 ## Here the resolution of the data is 'hours', we want to aggeregate the data up to 72 hours
 ## specify the aggregation durations

 durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)

 ## get the aggrageted data for each of the
station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",
 durations = durations)

 ## get initial values

 initial_params = egpd_idf_init(station_data = station_data,
             durations = durations, fitting_method = "mle",
               declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
               auto_fit = T, nrmse_tol = 0.1,use_r_optim = T, nrsme_quantile = 0)
 ## fit the data driven IDF
 fitted_idf = fit_egpd_idf_data_driven(station_data = station_data, durations = durations,
             declustering_duration = c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72),
             fitting_method = 'mle',  initial_params = initial_params,
              optim_algo = "BFGS")

 #check 'optim' params, for convergencem etc
 fitted_idf$fitted_params

 #parameters of the IDF
 fitted_idf$fitted_params$par

 # fitted egpd parameters for the given durations
 kappa_fit =  fitted_idf$kappa_param
 sigma_fit = fitted_idf$scale_param
 xi_fit = fitted_idf$shape_param

 #compute nrmse to check quality of fit
 nrmse_d = compute_nrsme(station_data, c(1,2,3,6,10,12, 16,18, 24, 48, 72),
  kappa_fit, sigma_fit, xi_fit, init_time_step = 1, q = 0)
 nrmse_d

# Plot the IDF curves
plot_egpdidf_curves(station_data = station_data,  kappa_fit = kappa_fit,
sigma_fit = sigma_fit, xi_fit = xi_fit, durations,
 declustering_duration=c(1,2,3,6,10,12, 16,18, 24, 48, 72), npy = 92, init_time_step=1 )
 
```

## Fitting IDF models that are based on scaling principle

This class of IDF models rely on the scaling of precipitation process in time. We consider eight different variations of this models.

The first four are all based on simple scaling model and some extensions

* Simple scaling IDF model
* Simple scaling IDF model with *scaling break*
* Simple scaling IDF model with *shape parameter as a function of duration*
* Simple scaling IDF model with *scaling break* and *shape parameter as a function of duration*

The last four are based on the General IDF formulation of Koutsoyiannis et al 1998.

* General IDF model
* General IDF model with *scaling break*
* General IDF model with *shape parameter as a function of duration*
* General IDF model with *scaling break* and *shape parameter as a function of duration*

For their details, refer to Haruna et al. 2020.

Here we give example of the simple-scaling IDF model. This is achieved by setting the  arguments  **simple_scaling = T, multi_regime = F, xi_constant = T**. The other seven models can be fitted by modifying these arguments. For example, in the case of IDF based on the General formulation of Koutsoyiannis 1998, the arguments should be **simple_scaling = F**, while **multi_regime = F, xi_constant = T**

```{r, eval=FALSE}

 ## load the data
 data("precipdata")

 ## Here the resolution of the data is in 'hours', we want to aggeregate the data up to 72 hours
 ## specify the aggregation durations

 durations =  c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72)

 ## get the aggrageted data for each of the
 station_data= aggregate_data(sample_data = precipdata, st_code = "SCH",  durations = durations)

 ## get initial values

 initial_params = egpd_idf_init(station_data = station_data,
                    durations = durations, fitting_method = "mle",
                   declustering_duration =  c(1,2,3,6,10,12, 16,18, 24, 48, 72),
                  auto_fit = T, nrmse_tol = 0.1,use_r_optim = T, nrsme_quantile = 0)
 
 ## fit the simple scaling  IDF model
 fitted_idf = fit_egpd_idf_scaling_models(station_data = station_data, durations = durations,
                 censored = initial_params$fits$lower_threshold,
                 declustering_duration = c(1,2, 3,  6,  10, 12,  16, 18,  24, 48, 72),
                 fitting_method = 'mle',  initial_params = initial_params,
              simple_scaling = T, multi_regime = F, xi_constant = T,
              init_time_step = 1,auto_fit = F, nrmse_tol = 0.1, use_r_optim = F,
              nrsme_quantile = 0, optim_algo = "BFGS")

 #check 'optim' params, for convergence etc
 fitted_idf$fitted_params

 #parameters of the IDF
 fitted_idf$fitted_params$par

 # fitted egpd parameters for the given durations
 kappa_fit =  fitted_idf$kappa_param
 sigma_fit = fitted_idf$scale_param
 xi_fit = fitted_idf$shape_param

 #compute nrmse to check quality of fit
 nrmse_d = compute_nrsme(station_data, c(1,2,3,6,10,12, 16,18, 24, 48, 72),
 kappa_fit, sigma_fit, xi_fit, init_time_step = 1, q = 0)
 nrmse_d

# Plot the IDF curves
plot_egpdidf_curves(station_data = station_data,  kappa_fit = kappa_fit,
 sigma_fit = sigma_fit, xi_fit = xi_fit, durations,
 declustering_duration=c(1,2,3,6,10,12, 16,18, 24, 48, 72), npy = 92, init_time_step=1 )
 

```


# References

- Naveau, P., Huser, R., Ribereau, P., and Hannart, A.: Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection,Water Resour. Res., 52, 2753–2769, https://doi.org/10.1002/2015WR018552, 2016.
- Haruna, Abubakar, Juliette Blanchet, and Anne-Catherine Favre. "Modeling Intensity-Duration-Frequency curves for the whole range of precipitation: A comparison of models." Authorea Preprints (2022)
- Koutsoyiannis, Demetris, Demosthenes Kozonis, and Alexandros Manetas. "A mathematical framework for studying rainfall intensity-duration-frequency relationships." Journal of hydrology 206.1-2 (1998): 118-135.
