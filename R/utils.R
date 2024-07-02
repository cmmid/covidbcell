#' A collection of utility functions for the analysis of the vaccine data
#' 
#' A nice theme for ggplot2 plots
david_theme_A <- function() {
    theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 14),
        title = element_text(size = 14),
       # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
    )
}

#' A nice theme for ggplot2 plots
david_theme_B <- function() {
    theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16),
        title = element_text(size = 16),
       # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
    )
}

#' A function to simulate the in-host kinetic model outside of stan
#' @param parms A list of parameters for the model
#' @param end_time The time to simulate the model to
#' @param init_cond The initial conditions for the model
#' @return A data frame of the model output
#' @export
ode_results_vac <- function(parms, end_time, init_cond) {

  des = function(t, state, parms) {
    with(as.list(c(state, parms)), {

     #   dB1 <- exp(- deltaAg_i * t) * param1 - B1 * param2a - B1 * deltaM_i
     #   dB2 <- B1 * param2a - exp(- deltaAg_i * t) * B2 * param2b - B2 * deltaM_i

        dB <- exp(- deltaAg * t) * param1 - exp(- deltaAg * t) * B * param2 - B * deltaM
        dP <- exp(- deltaAg * t) * B * param2 * prop_asc - P * deltaP
        dA1 <- P * param3 - A1 * deltaA
        dG <- exp(- deltaAg * t) * B * param2 * (1 - prop_asc) - G * deltaG
        dL <- G * deltaG - L * deltaL
        dA2 <- L * param3 * param4 - A2 * deltaA

        return(list(c(dB, dP, dA1, dG, dL, dA2)))
    })
  }

  yini = init_cond
  times = seq(0, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms)
  return(out)

}

#' A function to extract the results of the in-host kinetic model
#' @param params A list of parameters for the model
#' @param times The times to extract the results at
#' @param outputs The outputs to extract
#' @param init The initial conditions for the model
#' @return A data frame of the model output
#' @export
get_results_vac <- function(params, times, outputs, init) {
  t_max <- max(times)
  all_res <- ode_results_vac(params, t_max, init)
  actual_res <- all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[,outputs])
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = "")))
}

#' A function to calculate the combined mean of a set of means
#' @param means A vector of means
#' @param ns A vector of sample sizes
#' @return The combined mean
#' @export
calculate_combined_mean <- function(means, ns) {
  combined_mean <- sum(means * ns) / sum(ns)
  return(combined_mean)
}

#' A function to calculate the combined variance of a set of means
#' @param means A vector of means
#' @param ses A vector of standard errors
#' @param ns A vector of sample sizes
#' @return The combined variance
#' @export
calculate_combined_variance <- function(means, ses, ns) {
  variances <- ses^2 * ns
  n_total <- sum(ns)
  weighted_mean_diff <- sum(ns * (means - calculate_combined_mean(means, ns))^2)
  combined_variance <- (sum((ns - 1) * variances) + weighted_mean_diff) / (n_total - 1)
  return(combined_variance)
}

#' A function to convert confidence intervals to standard errors
#' @param ci_lowers A vector of lower confidence intervals
#' @param ci_uppers A vector of upper confidence intervals
#' @return A vector of standard errors
#' 
ci_to_se <- function(ci_lowers, ci_uppers) {
  ses <- (ci_uppers - ci_lowers) / (2 * 1.96)
  return(ses)
}

#' A function to calculate the combined lower confidence interval of a set of means
#' @param means A vector of means
#' @param ci_lowers A vector of lower confidence intervals
#' @param ci_uppers A vector of upper confidence intervals
#' @param ns A vector of sample sizes
#' @return The combined lower confidence interval
#' 
calculate_combined_ci_lower <- function(means, ci_lowers, ci_uppers, ns) {
  ses <- ci_to_se(ci_lowers, ci_uppers)
  combined_mean <- calculate_combined_mean(means, ns)
  combined_variance <- calculate_combined_variance(means, ses, ns)
  combined_se <- sqrt(combined_variance / sum(ns))
  ci_lower <- combined_mean - 1.96 * combined_se
  return(ci_lower)
}

#' A function to calculate the combined upper confidence interval of a set of means
#' @param means A vector of means
#' @param ci_lowers A vector of lower confidence intervals
#' @param ci_uppers A vector of upper confidence intervals
#' @param ns A vector of sample sizes
#' @return The combined upper confidence interval
#' 
calculate_combined_ci_upper <- function(means, ci_lowers, ci_uppers, ns) {
  ses <- ci_to_se(ci_lowers, ci_uppers)
  combined_mean <- calculate_combined_mean(means, ns)
  combined_variance <- calculate_combined_variance(means, ses, ns)
  combined_se <- sqrt(combined_variance / sum(ns))
  ci_upper <- combined_mean + 1.96 * combined_se
  return(ci_upper)
}

#' A function to exttact the data from a stan model
#' @param foldername The name of the folder containing the stan model
#' @return A data frame of the stan model output
#' @export
get_data_ind <- function(foldername) {
    data_model <- readRDS(here::here("outputs", "data_clean", foldername, "df_meta.RDS"))
    data_listA <- readRDS(here::here("outputs", "data_clean", foldername, "list_stan_data.RDS"))
    PIDs <- data_model$PID
    N <- length(PIDs)

    data_i_all <- 
        map_df( 1:N,
        function(ind) {
            data_listA$y_obs[[ind]] %>% as.data.frame %>% setNames(c("B", "P", "A")) %>%
                mutate(t = data_listA$y_times[ind, ]) %>%
                pivot_longer(!t, names_to = "state", values_to = "value") %>% filter(value > -1) %>% mutate(PID = PIDs[ind])
        }
    )
    relabel_state_variable <- c("A" = "sVNT to ancestral variant", "P" = "Plasmablast conc.", "B" = "Memory B-cell conc.")
    data_i_all_labs <- data_i_all %>% mutate(state = recode(state, !!!relabel_state_variable) )

}

# A function to de-anonymize the PID
#' @param df A data frame containing the PID
#' @return A data frame with the PID de-anonymized
#' @export
deanonymize_pid <- function(df) {
    df %>% mutate(PID = factor(PID, levels = unique(PID))) %>%
        mutate(PID = as.numeric(PID)) %>% 
        mutate(PID = factor(PID, levels = as.character(1:41)))
}

#' Some cute colours
color_study_A <- "#fd766a"
color_study_B <- "#5b84b1"


color_MBC <- "#89c462"
color_PB <- "#1cae73"
color_sVNT <- "#18abb0"