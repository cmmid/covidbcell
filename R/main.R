# Load this file first to get all the necessary packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, ggplot2, here, posterior, tidybayes,
    cmdstanr, rstan, deSolve, patchwork, Rcpp, devtools, ptmc, coda, bayesplot,
    extraDistr, hmer, loo, ggh4x)

#source(here::here("R/utils.R"))

inv.logit <- function(x) {
    exp(x)/(1+exp(x)) 
} 

# Load all the functions

# Load the helper functions
source("R/utils.R")
# Load the data cleaning functions
source("R/clean_data.r")
# Load the data cleaning plot functions
source("R/clean_data_plot.r")
# Load the function to assess the quality of the mcmc fit
source("R/met_mcmc.R")
# Load the prior predictive check functions
source("R/met_pp_check.R")
# Load the functions to plot the posterior predictive distributions of the antibody kinetics
source("R/res_abkin.R")
# Load the functions to plot the posterior predictive distributions of the calibration data set fit
source("R/res_fit_cal.R")
# Load the functions to plot the posterior predictive distributions of the validation data set fit
source("R/res_fit_val.R")
# Load the functions to plot the raw posteriors of the fitted model
source("R/res_post_desc.R")

