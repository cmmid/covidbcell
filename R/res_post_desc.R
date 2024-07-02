
#' Plot the raw posterior distributions for each of the parameters
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
plot_posteriors <- function(stanfit, fig_folder) {    
    
    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike model"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD model"
        antigen_file <- "rbd"
    }

    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))

    recode_expo <- c("2" = "ChAdOx1", "1" = "BNT162b2")
    recode_time <- c("1" = "<=28 days", "2" = "28+ days")
    recode_age <- c("1" = "<30 years", "2" = "30–39 years", "3" = "40–49 years", "4" = "50–59 years", "5" = "60+ years")

    param1_post <- bind_rows(
        stanfit$draws(c("param1", "v_1", "sigmav1")) %>% spread_draws(param1, v_1[x], sigmav1) %>% 
            mutate("theta_1" = v_1 * sigmav1, .keep = "unused") %>% 
            mutate(covar = "Vaccine type") %>% mutate(x = recode(x, !!!recode_expo)),
        stanfit$draws(c("param1", "t_1", "sigmat1")) %>% spread_draws(param1, t_1[x], sigmat1) %>% 
            mutate("theta_1" = t_1 * sigmat1, .keep = "unused") %>% mutate(covar = "Time since first dose") %>%
            mutate(x = recode(x, !!!recode_time)),
        stanfit$draws(c("param1", "a_1", "sigmaa1")) %>% spread_draws(param1, a_1[x], sigmaa1) %>% 
            mutate("theta_1" =  a_1 * sigmaa1, .keep = "unused") %>% mutate(covar = "Age group") %>%
            mutate(x = recode(x, !!!recode_age))
    ) %>% mutate(x = factor(x, levels = c(recode_age, recode_time, recode_expo)))
    

    param3_post <- bind_rows(
        stanfit$draws(c("param3", "v_3", "sigmav3")) %>% spread_draws(param3, v_3[x], sigmav3) %>% 
            mutate("theta_3" = v_3 * sigmav3, .keep = "unused") %>% 
            mutate(covar = "Vaccine type") %>% mutate(x = recode(x, !!!recode_expo)),
        stanfit$draws(c("param3", "t_3", "sigmat3")) %>% spread_draws(param3, t_3[x], sigmat3) %>% 
            mutate("theta_3" = t_3 * sigmat3, .keep = "unused") %>% mutate(covar = "Time since first dose") %>%
            mutate(x = recode(x, !!!recode_time)),
        stanfit$draws(c("param3", "a_3", "sigmaa3")) %>% spread_draws(param3, a_3[x], sigmaa3) %>% 
            mutate("theta_3" =  a_3 * sigmaa3, .keep = "unused") %>% mutate(covar = "Age group") %>%
            mutate(x = recode(x, !!!recode_age))
    ) %>% mutate(x = factor(x, levels = c(recode_age, recode_time, recode_expo)))

    param4_post <- bind_rows(
        stanfit$draws(c("param4", "v_4", "sigmav4")) %>% spread_draws(param4, v_4[x], sigmav4) %>% 
            mutate("theta_4" = v_4 * sigmav4, .keep = "unused") %>% 
            mutate(covar = "Vaccine type") %>% mutate(x = recode(x, !!!recode_expo)),
        stanfit$draws(c("param4", "t_4", "sigmat4")) %>% spread_draws(param4, t_4[x], sigmat4) %>% 
            mutate("theta_4" = t_4 * sigmat4, .keep = "unused") %>% mutate(covar = "Time since first dose") %>%
            mutate(x = recode(x, !!!recode_time)),
        stanfit$draws(c("param4", "a_4", "sigmaa4")) %>% spread_draws(param4, a_4[x], sigmaa4) %>% 
            mutate("theta_4" =  a_4 * sigmaa4, .keep = "unused") %>% mutate(covar = "Age group") %>%
            mutate(x = recode(x, !!!recode_age))
    ) %>% mutate(x = factor(x, levels = c(recode_age, recode_time, recode_expo)))


    p1 <- param1_post %>%
            ggplot() + 
                stat_pointinterval(aes(x = covar, y = theta_1, color = x), position = position_dodge(0.5)) + 
                theme_bw() + labs(x = "Covariate", y = "Posterior distribution", color = "Covariate levels", 
                title = "Rate of B-cell proliferation per vaccine dose")
    p3 <- param3_post %>%
            ggplot() + 
                stat_pointinterval(aes(x = covar, y = theta_3, color = x), position = position_dodge(0.5)) + 
                theme_bw() + labs(x = "Covariate", y = "Posterior distribution", color = "Covariate levels",
                title = "Rate of antibody production per plasmablast")
    p4 <- param4_post %>%
            ggplot() + 
                stat_pointinterval(aes(x = covar, y = theta_4, color = x), position = position_dodge(0.5)) + 
                theme_bw() + labs(x = "Covariate", y = "Posterior distribution", color = "Covariate levels",
                title = "Rate of antibody production per plasma cell")
    p1 / p3 / p4 + plot_layout(guides = "collect") + plot_annotation(title = antigen_str,  theme = theme(
      plot.title = element_text(size = 24) )) 
    ggsave(here::here("outputs", "figs", "post_desc", paste0("post_only_" , antigen_file, ".png")), height = 10, width = 12)

}

#' Plot the posterior distributions for the decay parameters of the model
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' @return A ggplot object
plot_posteriors_decay <- function(stanfit, fig_folder) {    

    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike model"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD model"
        antigen_file <- "rbd"
    }

    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))

    recode_expo <- c("2" = "ChAdOx1", "1" = "BNT162b2")
    recode_time <- c("1" = "<=28 days", "2" = "28+ days")
    recode_age <- c("1" = "<30 years", "2" = "30–39 years", "3" = "40–49 years", "4" = "50–59 years", "5" = "60+ years")

    require(tidybayes)
    decays_names <- c("deltaAg", "deltaP", "deltaL", "deltaA", "deltaG")
    df_post_decay <- stanfit$draws(decays_names) %>% as_draws_df %>% mutate(type = "Posterior distribution") %>% mutate(deltaAg = deltaAg * 30)
    source(here::here("R", "met_pp_check.R"))
    df_prior_decay <- create_priors_A() %>% select(decays_names) %>% mutate(type = "Prior distribution") %>%
        mutate(deltaAg = 1 / (deltaAg) , deltaP = 1 / deltaP, deltaL = 1 / deltaL, deltaA = 1 / deltaA, deltaG = 1 / deltaG)

    df_compare_decay <- bind_rows(
        df_post_decay,
        df_prior_decay
    ) %>%
        pivot_longer(decays_names)

    df_compare_decay %>% ggplot() + 
        geom_violin(aes(y = type, x = value)) + 
        facet_wrap(vars(name), scales = "free") + theme_bw() + 
        labs(x = "Mean rate of exponential decay for each biomarker", y = "") + plot_annotation(title = antigen_str,  theme = theme(
      plot.title = element_text(size = 24) )) 
    ggsave(here::here("outputs", "figs", "post_desc", paste0("post_decay_" , antigen_file, ".png")), height = 10, width = 12)

}

#' Plot the marginal posterior distributions for the hierarchicial effects in the model 
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
plot_marg_posteriors <- function(stanfit, fig_folder) {

    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike model"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD model"
        antigen_file <- "rbd"
    }

    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))

    #pp_ll <- post_sample_cross$draws(c("lp__"))  %>% ggplot() + geom_line(aes(sample_no, lpost, color = chain_no))
    #ggsave(here::here("outputs", "figs", "pp_ll.pdf"), height = 20, width = 20)

    ################################################
    ## Plot postiors of the regression effects ##
    ################################################

    recode_expo <- c("1" = "ChAdOx1", "2" = "BNT162b2")
    recode_time <- c("1" = "<=28 days", "2" = "28+ days")
    recode_age <- c("1" = "<30 years", "2" = "30–39 years", "3" = "40–49 years", "4" = "50–59 years", "5" = "60+ years")
    recode_states <- c("v" = "Exposure type", "t" = "Time since first dose", "a" = "Age group at second dose")
    relabel_model_state <- c("1" = "a1: Rate of \nB-cell proliferation",
        "2" = "a2: Rate of derivation \nto plasmablasts", 
        "3" = "a3: Rate of antibody \nproduction per plasmablast", 
        "4" = "a4: Rate of antibody \nproduction per plasma cell")

    p1 <- stanfit$draws(c("theta_v")) %>% spread_draws(theta_v[i, s]) %>% 
        pivot_wider(names_from = "s", values_from = "theta_v") %>%
        mutate(`4` = `3` * `4`) %>%
        pivot_longer(`1`:`4`, names_to = "s", values_to = "theta_v") %>%
        mutate(s = recode(s, !!!relabel_model_state)) %>%
        mutate(i = recode(i, !!!recode_expo)) %>%
        ggplot() + stat_pointinterval(aes(x = s, y = theta_v,
            color = as.character(i)), position = position_dodge(0.5)) + 
            labs(x = "Parameter value", y = "Rate of B-cell proliferation per dose", color = "Vaccine type")
            
    p2 <- stanfit$draws(c("theta_t")) %>% spread_draws(theta_t[i, v, s]) %>% 
        pivot_wider(names_from = "s", values_from = "theta_t") %>%
        mutate(`4` = `3` * `4`) %>%
        pivot_longer(`1`:`4`, names_to = "s", values_to = "theta_t") %>%
        mutate(v = recode(v, !!!recode_expo)) %>%
        mutate(s = recode(s, !!!relabel_model_state)) %>%
        mutate(i = recode(i, !!!recode_time)) %>%
        ggplot() + stat_pointinterval(aes(x = s, y = theta_t,
            color = as.character(i)), position = position_dodge(0.5)) + 
            labs(x = "Parameter value", y = "Rate of B-cell proliferation per dose", color = "Time since first dose") + 
            facet_grid(vars(v))

    p3 <- stanfit$draws(c("theta_a")) %>% spread_draws(theta_a[i, v, s]) %>% 
        pivot_wider(names_from = "s", values_from = "theta_a") %>%
        mutate(`4` = `3` * `4`) %>%
        pivot_longer(`1`:`4`, names_to = "s", values_to = "theta_a") %>%
        mutate(v = recode(v, !!!recode_expo)) %>%
        mutate(s = recode(s, !!!relabel_model_state)) %>%
        mutate(i = recode(i, !!!recode_age)) %>%
        ggplot() + stat_pointinterval(aes(x = s, y = theta_a,
            color = as.character(i)), position = position_dodge(0.5)) + 
            labs(x = "Parameter value", y = "Rate of B-cell proliferation per dose", color = "Age group") + 
            facet_grid(vars(v))

    p1 / p2 / p3 + plot_layout(guides = "collect") + plot_annotation(title = antigen_str,  theme = theme(
      plot.title = element_text(size = 24) )) & theme_bw()
    ggsave(here::here("outputs", "figs", "post_desc", paste0("marg_post_only_" , antigen_file, ".png")), height = 10, width = 12)

}

#' Plot the posterior distributions for the B cell kinetics model
#' @param stanfit The stanfit object
#' @param fig_folder The folder where the data is stored
#' @return A ggplot object
plot_bcell <- function(stanfit, fig_folder) {

    if (fig_folder == "nih_wu_s") {
        antigen_str <- "Ancestral spike model"
        antigen_file <- "s"
    } else {
        antigen_str <- "Ancestral RBD model"
        antigen_file <- "rbd"
    }

    data_model <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS"))

    recode_expo <- c("2" = "ChAdOx1", "1" = "BNT162b2")

    relabel_model_state <- c("1" = "a1: Rate of B-cell proliferation",
        "2" = "a2: Rate of derivation to plasmablasts", 
        "3" = "a3: Rate of neutralising antibody production from plasmablasts", 
        "4" = "a4: Rate of neutralising antibody production from plasma cells")

    comparison_t <- stanfit$draws(c("theta_t")) %>%
        spread_draws(theta_t[i, v, s]) %>% 
            mutate(s = recode(s, !!!relabel_model_state)) %>%
            mutate(v = recode(v, !!!recode_expo)) %>%
            filter(s != "a2: Rate of derivation to plasmablasts") %>% 
            pivot_wider(names_from = s, values_from = theta_t) 
    comparison_t_mean <- comparison_t %>% summarise(
            `a1: Rate of B-cell proliferation` = median(`a1: Rate of B-cell proliferation`),
            `a3: Rate of neutralising antibody production from plasmablasts` = median(`a3: Rate of neutralising antibody production from plasmablasts`),
            `a4: Rate of neutralising antibody production from plasma cells` = median(`a4: Rate of neutralising antibody production from plasma cells`)) %>%
            mutate(i = recode(i, "1" = "<28 days", "2" = "\u2265 28 days")) 

    p1 <- comparison_t %>%  
            ggplot() + 
            geom_point(aes(x = `a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, color = v
            ), 
                size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_t_mean, aes(`a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_t_mean, aes(`a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Time since first dose")

    p2 <- comparison_t %>%  
            ggplot() + 
            geom_point(aes(x = `a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, color = v
            ),  size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_t_mean, aes(`a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_t_mean, aes(`a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Time since first dose")

    p3 <- comparison_t %>%
            ggplot() + 
            geom_point(aes(x = `a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, color = v
            ),  size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_t_mean, aes(`a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_t_mean, aes(`a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Time since first dose")

    p_time_spike <- wrap_elements(((p1 | p2 | p3) + plot_layout(guides = "collect") + plot_annotation(
        title = paste0("Impact of time since first dose on humoral kinetics of ", antigen_str))) & guides(size = "none", color = "none", 
                fill = guide_legend(override.aes = list(shape=16, size = 5, color = c("#2271B2", "#D55E00")))
                ) )

    saveRDS(p_time_spike, here::here("outputs", "figs", "bcell", paste0("time_" , antigen_file, ".RDS")))
    ggsave(here::here("outputs", "figs", "bcell", paste0("time_" , antigen_file, ".png")), height = 7, width = 15)


    comparison_a <- stanfit$draws(c("theta_a")) %>%
        spread_draws(theta_a[i, v, s]) %>% 
            mutate(s = recode(s, !!!relabel_model_state)) %>%
            mutate(v = recode(v, !!!recode_expo)) %>%
            filter(s != "a2: Rate of derivation to plasmablasts") %>% 
            pivot_wider(names_from = s, values_from = theta_a) 
    comparison_a_mean <- comparison_a %>% summarise(
            `a1: Rate of B-cell proliferation` = median(`a1: Rate of B-cell proliferation`),
            `a3: Rate of neutralising antibody production from plasmablasts` = median(`a3: Rate of neutralising antibody production from plasmablasts`),
            `a4: Rate of neutralising antibody production from plasma cells` = median(`a4: Rate of neutralising antibody production from plasma cells`)) %>% 
            mutate(i = recode(i, "1" = "<30", "2" = "30–39", "3" = "40–49", "4" = "50–59", "5" = "60+")) 

    p1 <- comparison_a %>%  
            ggplot() + 
            geom_point(aes(x = `a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, color = v
            ), 
                size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_a_mean, aes(`a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_a_mean, aes(`a1: Rate of B-cell proliferation`, `a3: Rate of neutralising antibody production from plasmablasts`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Age groups (yrs)")


    p2 <- comparison_a %>%  
            ggplot() + 
            geom_point(aes(x = `a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, color = v
            ),  size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_a_mean, aes(`a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_a_mean, aes(`a1: Rate of B-cell proliferation`, `a4: Rate of neutralising antibody production from plasma cells`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Age groups (yrs)")

    p3 <- comparison_a %>%
            ggplot() + 
            geom_point(aes(x = `a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, color = v
            ),  size = 0.5, alpha = 0.2) + 
            geom_point(data = comparison_a_mean, aes(`a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, 
                shape = i), 
                size = 8, fill = "white", color = "white") +
            geom_point(data = comparison_a_mean, aes(`a3: Rate of neutralising antibody production from plasmablasts`, `a4: Rate of neutralising antibody production from plasma cells`, 
                fill = v, shape = i), 
                size = 6, alpha = 0.8) +
            scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
            scale_color_manual(values = c("#2271B2", "#D55E00")) +
            scale_fill_manual(values = c("#2271B2", "#D55E00")) +
            theme_bw() + labs(fill = "Vaccine type", shape = "Age groups (yrs)")


    p_age_spike <- wrap_elements(((p1 | p2 | p3) + plot_layout(guides = "collect") + plot_annotation(
        title = paste0("Impact of age on humoral kinetics of ", antigen_str))) & guides(size = "none", color = "none", 
                fill = guide_legend(override.aes = list(shape=16, size = 5, color = c("#2271B2", "#D55E00")))
                )  )
    saveRDS(p_age_spike, here::here("outputs", "figs", "bcell", paste0("age_" , antigen_file, ".RDS")))
    ggsave(here::here("outputs", "figs", "bcell", paste0("age_" , antigen_file, ".png")), height = 7, width = 15)

}