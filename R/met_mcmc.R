# Purpose: Functions for MCMC diagnostics and convergence checks
#' @title Make the plots showing the PRSF and ESS for the antigen-specific parameters
#' @param fit_model The stanfit object
#' @param title_str The title of the plot
#' @param param_vec The parameters to plot
#' @param font_size The font size for the x-axis
#' @return A ggplot object
#' @export
#' @examples
#' get_cc(fit_model, "Dynamic model parameters", fitparlist, 10)
#' get_cc(fit_model,"Hierarchical model parameters on a1, host factors", fithierlist_a1, 10)
#' get_cc(fit_model,"Hierarchical model parameters on a1, individual-level factors", fithierlist_a1_ind, 5)
#'  
get_cc <- function(fit_model, title_str, param_vec, font_size) {
    draws <- as_draws_df(fit_model$draws(param_vec))
    paras_vec_names <- names(draws)

    df_cc <- fit_model$summary(variables = param_vec) %>%
        select(variable, PRSF = rhat, ESS_bulk = ess_bulk, ESS_tail = ess_tail) %>%
            pivot_longer(!variable, names_to = "stat", values_to = "value") %>%
            mutate(stat = factor(stat, levels = c("PRSF", "ESS_bulk", "ESS_tail"))) %>%
            mutate(variable = factor(variable, levels = paras_vec_names))

    thresholds <-
        data.frame(
            PRSF = 1.1,
            ESS_bulk = 2000,
            ESS_tail = 2000
        ) %>% pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
        mutate(stat = factor(stat, levels = c("PRSF", "ESS_bulk", "ESS_tail")))


    df_cc %>%
        ggplot() +
            geom_hline(data = thresholds, aes(yintercept = value), linetype = "dashed", color = "gray") +
            geom_point(aes(variable, value), size = 4, shape = 21, fill = "gray") +
            facet_wrap(vars(stat), scales = "free") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, size = font_size)) +
            david_theme_A() +
            ggtitle(paste0(title_str))
}

#' @title Make all the plots showing the PRSF and ESS for the antigen-specific parameters
#' @param stanfit_name The name of the stanfit object
#' @param str_antigen The antigen to plot
#' @return A ggplot object
#' @export
run_cc_antigen <- function(stanfit_name, str_antigen) {

    fitparlist <- c("param1" , "param2", "param3", "param4", "prop_asc", "deltaAg", "deltaP", "deltaA", "deltaG", "deltaL", "sigma_M", "sigma_P", "sigma_A")
    fithierlist_a1 <- c("v_1" , "t_1", "a_1", "sigmav1", "sigmat1", "sigmaa1")
    fithierlist_a1_ind <- c("ind_1" , "sigmai1")

    fithierlist_a3 <- c("v_3" , "t_3", "a_3", "sigmav3", "sigmat3", "sigmaa3")
    fithierlist_a3_ind <- c("ind_3" , "sigmai3")

    fithierlist_a4 <- c("v_4" , "t_4", "a_4", "sigmav4", "sigmat4", "sigmaa4")
    fithierlist_a4_ind <- c("ind_3" , "sigmai3")


    fit_model <- readRDS(file = here::here("outputs", "stanfit",  paste0(stanfit_name, ".RDS")) )
    data_model <- readRDS(here::here("outputs", "data_clean", "nih_wu_s", "df_meta.RDS"))

    p0 <- get_cc(fit_model, "Dynamic model parameters", fitparlist, 10)

    p1 <- get_cc(fit_model,"Hierarchical model parameters on a1, host factors", fithierlist_a1, 10)
    p2 <- get_cc(fit_model,"Hierarchical model parameters on a1, individual-level factors", fithierlist_a1_ind, 5)
    p1A <- p1 / p2

    p1 <- get_cc(fit_model,"Hierarchical model parameters on a3, host factors", fithierlist_a3, 10)
    p2 <- get_cc(fit_model,"Hierarchical model parameters on a3, individual-level factors", fithierlist_a3_ind, 5)
    p1B <- p1 / p2


    p1 <- get_cc(fit_model,"Hierarchical model parameters on a4, host factors", fithierlist_a4, 10)
    p2 <- get_cc(fit_model,"Hierarchical model parameters on a4, individual-level factors", fithierlist_a4_ind, 5)
    p1C <- p1 / p2

    if (str_antigen == "s" ) {
        str_title <- "spike"
    } else if (str_antigen == "rbd") {
        str_title <- "RBD"
    }

    p0 / p1A / p1B / p1C + plot_layout(heights = c(1, 4, 4, 4)) + plot_annotation(title = paste0("Chain convergence and resolution for Ancestral ", str_title, "model"), 
        theme = theme(
      plot.title = element_text(size = 24) ))
    ggsave(here::here("outputs", "figs", "cc", paste0("cc_", str_antigen, ".png")), width = 12, height = 20)

}

#' @title Function to print the sampler diagnostics
#' @param fit_model The stanfit object
#' export
run_sampler_diag <- function(fit_model) {
    sampler_diagnostics <- fit_model$sampler_diagnostics()
    df_sampler_diagnostics <- as_draws_df(sampler_diagnostics) 

    paste0("There are ", df_sampler_diagnostics$divergent__ %>% sum, " divergent transitions. There are ", (df_sampler_diagnostics$treedepth__ > 10) %>% sum, ", which hit max tree depth of 10.")
}

#' @title Function to plot the trace and pairs plots to determine divergence transitions
#' @param fit_model The stanfit object
#' @param param_vec The parameters to plot
#' @param name_str The name of the plot
#' @export
#' 
plot_diag <- function(fit_model, param_vec, name_str) {

    draws <- as_draws_df(fit_model$draws(param_vec))
    param_vec_names <- setdiff(names(draws), c(".chain", ".iteration", ".draw"))

    p1 <- mcmc_pairs(
        draws,
        np = nuts_params(fit_model),
        off_diag_args = list(size = 1, alpha = 1/3),
        condition = pairs_condition(nuts = "accept_stat__"),
        pars = param_vec_names  # Replace with your parameter names
    )
    ggsave(here::here("outputs", "figs", "diag", paste0(name_str, "_pairs.png")), p1, width = 10, height = 10)
    p2 <- mcmc_trace(
        draws,
        np = nuts_params(fit_model),
        pars = param_vec_names  # Replace with your parameter names
    )
    ggsave(here::here("outputs", "figs", "diag", paste0(name_str, "_hist.png")), p2, width = 10, height = 10)
}
