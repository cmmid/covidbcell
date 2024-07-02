#' @title: Prior predictive checks
#' @description: Functions to generate prior predictive checks
#' This is the raw prior predictive checks without the logistic trnasformed parameters
create_priors_A <- function() {
    a1 <- runif(1000, 0, 2)
    a2 <- runif(1000, 0, 1)
    a3 <- runif(1000, 0, 6)
    a4 <- runif(1000, 0, 0.5)

    prop_asc <- rbeta(1000, 60, 40)

    deltaAg <- runif(1000, 1 / 30, 1)
    deltaP <- rnorm(1000, 4, 1)
    deltaA <- rnorm(1000, 30, 5)
    deltaG <- rnorm(1000, 18, 2)
    deltaL <- rnorm(1000, 730, 200)
    deltaM <- 0.001

    data.frame(
        param1 = a1,
        param2 = a2,
        param3 = a3,
        param4 = a4,
        prop_asc = prop_asc,
        deltaAg = 1 / (deltaAg * 30),        
        deltaM = deltaM,
        deltaP = 1 / deltaP,
        deltaL = 1 / deltaL,
        deltaA = 1 / deltaA,
        deltaG = 1/ deltaG
    )
}

#' @title: Prior predictive checks
#' @description: Functions to generate prior predictive checks
#' This is the prior predictive checks including the logistic trnasformed parameters
create_priors_B <- function(exp_val) {

    b1 <- rnorm(1000, 0, 1.81)
    b2 <- rnorm(1000, 0, 1.81)
    b3 <- rnorm(1000, 0, 1.81)
    b4 <- rnorm(1000, 0, 1.81)

    # hierarchcial effects
    noncen <- function() {
        z <- rnorm(1000, 0, 1)
        sigma <- rexp(1000, exp_val)
        z * sigma
    }

    a1_t <- b1 + noncen() + noncen() + noncen() + noncen();
    a2_t <- b2;
    a3_t <- b3 + noncen() + noncen() + noncen() + noncen();
    a4_t <- b4 + noncen() + noncen() + noncen() + noncen();

    a1 <- inv.logit(a1_t) * 2
    a2 <- inv.logit(a2_t)
    a3 <- inv.logit(a3_t) * 6
    a4 <- inv.logit(a4_t) * 0.5

    prop_asc <- rbeta(1000, 60, 40)

    deltaAg <- runif(1000, 1/30, 1)
    deltaP <- rnorm(1000, 4, 1)
    deltaA <- rnorm(1000, 30, 5)
    deltaG <- rnorm(1000, 18, 2)
    deltaL <- rnorm(1000, 730, 200)
    deltaM <- 0.001

    data.frame(
        param1 = a1,
        param2 = a2,
        param3 = a3,
        param4 = a4,
        prop_asc = prop_asc,
        deltaAg = 1 / (deltaAg * 30),
        deltaM = deltaM,
        deltaP = 1 / deltaP,
        deltaL = 1 / deltaL,
        deltaA = 1 /  deltaA,
        deltaG = 1/ deltaG
    )
}

#' Function to find the c value for the prior predictive checks
create_find_cvalue <- function() {
    prior_model_A <- create_priors_A() %>% mutate(model = "Model A" )
    prior_model_B_1 <- create_priors_B(1) %>% mutate(model = "Model B" )
    prior_model_B_2 <- create_priors_B(3) %>% mutate(model = "Model B")

     bind_rows(prior_model_A, prior_model_B_1) %>% 
        select(param1, param2, param3, param4, model) %>% pivot_longer(!model) %>%
        mutate(name = recode(name, param1 = "a_1", param2 = "a_2", param3 = "a_3", param4 = "a_4")) 


    p1 <- bind_rows(prior_model_A, prior_model_B_1) %>% 
        select(param1, param2, param3, param4, model) %>% pivot_longer(!model) %>%
        mutate(name = recode(name, param1 = "a_1", param2 = "a_2", param3 = "a_3", param4 = "a_4")) %>%
        ggplot() + geom_violin(aes(x = name, y = value, fill = model), position = position_dodge(0.8)) +
        ggtitle("Prior predictive distributions for each  when pooling prior is exponential(1)") + theme_bw() + 
        scale_fill_manual(
            values = c("#5b84b1", "#be1500"),
            labels = c(
            expression(U(0, b[i])),
            expression(g(theta[i]))
            )
        ) + 
        labs( x = "Parameter", y = "Value", fill = "Prior model")

    p2 <- bind_rows(prior_model_A, prior_model_B_2) %>% 
        select(param1, param2, param3, param4, model) %>% pivot_longer(!model) %>%
        mutate(name = recode(name, param1 = "a_1", param2 = "a_2", param3 = "a_3", param4 = "a_4")) %>%
                ggplot() + geom_violin(aes(x = name, y = value, fill = model), position = position_dodge(0.8)) +
        ggtitle("Prior predictive distributions for each  when pooling prior is exponential(3)")  + theme_bw() + 
                scale_fill_manual(
            values = c("#5b84b1", "#be1500"),
            labels = c(
            expression(U(0, b[i])),
            expression(g(theta[i]) ))
            ) +
            labs( x = "Parameter", y = "Value", fill = "Prior model")



    p1 / p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position='top')
    ggsave(here::here("outputs", "figs", "prior", "find_c.png"), width = 10, height = 10)
}

#' Function to plot the prior predictive checks (figure SM5)
plot_prior_pred <- function() {

    data_meta <- readRDS(here::here("outputs", "data_clean", fig_folder, "df_meta.RDS")) %>% 
        mutate(time_d1_i = recode(time_d1, "<28 days" = "1", "28+ days" = "2")) %>% 
        mutate(vax_type_i = recode(vax_type, "BNT162b2" = "2", "ChAdOx1" = "1"))

    data_listA <- readRDS(here::here("outputs", "data_clean", fig_folder, "list_stan_data.RDS"))

    PIDs <- data_meta$PID

    prior_model <- create_priors_A()
    M <- 100
    df_uncert_all <- map_df(1:length(PIDs), 
        function(s) {
            init_i <- data_listA$y0_init[s, ]
            names(init_i) <- c("B", "P", "A1", "G", "L", "A2")

            df_uncert <- map_df(1:M,
                ~ode_results_vac(prior_model[.x, ], 220, init_i) %>% as.data.frame  )  %>%
                select(time, B, P, A1, A2) %>% mutate(A = A1 + A2) %>% 
                mutate(i = map(1:M, ~rep(.x, 221)) %>% unlist) %>% 
                mutate(PID = PIDs[s])
        }
    )

    df_uncert_all_clean <- df_uncert_all %>% 
        pivot_longer(c(B, P, A), names_to = "biomarker", values_to = "value") %>% 
        mutate(biomarker = recode(biomarker, B = "Memory B cells", P = "Plasmablasts", A = "sVNT")) %>% 
        mutate(biomarker = factor(biomarker, levels = c("Memory B cells", "Plasmablasts", "sVNT"))) %>% 
        mutate(PID = factor(PID, levels = PIDs))

    df_uncert_all_clean_sum <- df_uncert_all_clean %>% group_by(biomarker, time) %>% mean_qi(value) 


    df_data_s <- readRDS(here::here("outputs", "data_clean", "nih_wu_s", "df_model_data.RDS"))
    df_data_rbd <- readRDS(here::here("outputs", "data_clean", "nih_wu_rbd", "df_model_data.RDS"))

    get_plot_data <- function(df_data) {
        df_data %>% 
        select(PID, time_post_exp, Bmem, pbBmem, IC50hs) %>% 
        pivot_longer(c(Bmem, pbBmem, IC50hs), names_to = "biomarker", values_to = "value") %>%
        mutate(biomarker = recode(biomarker, Bmem = "Memory B cells", pbBmem = "Plasmablasts", IC50hs = "sVNT")) %>% 
        filter(value > -1)
    }

    df_data_clean <- bind_rows(
        get_plot_data(df_data_s) %>% mutate(antigen = "Ancestral spike"),
        get_plot_data(df_data_rbd)  %>% mutate(antigen = "Ancestral spikRBDe")
    )

    df_data_clean %>%
        ggplot() + 
            geom_point(aes(x = time_post_exp, y = value, shape = antigen, color = biomarker), alpha = 0.5) +
            facet_wrap(vars(biomarker), scales = "free_y") + 
            geom_ribbon(data = df_uncert_all_clean_sum, aes(x = time, ymin = .lower, ymax = .upper, fill = biomarker), alpha = 0.3) + 
            geom_line(data = df_uncert_all_clean_sum,aes(x = time, y = value, color = biomarker)) + 
            scale_color_manual(values = c(color_MBC, color_PB, color_sVNT)) +
            scale_fill_manual(values = c(color_MBC, color_PB, color_sVNT)) +
            theme_bw() + 
            theme(legend.position = "none") + 
            labs(x = "Time after second dose (days)", y = "Biomarker value") + 
            ggtitle("Prior predictive trajectories for each biomarker") + 
            david_theme_A() 
    ggsave(here::here("outputs", "figs", "prior", "prior_pred.png"), width = 10, height = 8)
}