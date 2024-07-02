#' A function to make figure SM2 in Supplementary Material
plot_data_clean <- function(data) {
    df_nih_s <- readRDS(here::here("outputs", "data_clean", "nih_wu_s", "df_model_data.RDS")) %>% mutate(data = "Study A (Calibration)")
    df_nih_rbd <- readRDS(here::here("outputs", "data_clean", "nih_wu_rbd", "df_model_data.RDS"))  %>% mutate(data = "Study A (Calibration)")
    df_sid_s <- readRDS(here::here("outputs", "data_clean", "sid_wu_s", "df_model_data.RDS"))  %>% mutate(data = "Study B (Validation)")
    df_sid_rbd <- readRDS(here::here("outputs", "data_clean", "sid_wu_rbd", "df_model_data.RDS"))  %>% mutate(data = "Study B (Validation)")

    data_compare <- bind_rows(df_nih_s, df_nih_rbd, df_sid_s, df_sid_rbd)  %>%
        pivot_longer(c(Bmem, pbBmem, IC50hs), names_to = "variable", values_to = "value") %>%
        filter(value > -1) %>% 
        mutate(biomarker = recode(variable, Bmem = "Memory B cell", pbBmem = "Plasmablasts", IC50hs = "sVNT")) %>%
        mutate(biomarker = factor(biomarker, levels = c("Memory B cell", "Plasmablasts", "sVNT")))

    # Plot of initial values

    p1 <- data_compare %>% filter(VisitN == 0) %>%
        ggplot() + 
            geom_point(aes(x = value, y = "", color = data, shape = data), alpha = 0.6, size = 3, position = position_dodge(0.8)) + 
            geom_boxplot(aes(x = value, y = "", color = data, shape = data), alpha = 0.6, size = 1, position = position_dodge(0.8)) + 
            scale_color_manual(values = c("#fd766a", "#5b84b1")) + 
            facet_grid(cols = vars(biomarker), rows = vars(antigen), scales = "free_x") + 
            theme_bw() + 
            theme( strip.background = element_blank(), legend.position = "top") +
            labs(x = "Biomarker value", y = "Antigen", color = "Dataset", shape = "Dataset") + 
            david_theme_A() + 
            ggtitle("Biomarker data priors to second dose (baseline)")

    p2 <- data_compare %>% filter(variable != "IC50hs") %>% 
        ggplot() +
        geom_smooth(aes(x = time_post_exp, y = value), se = FALSE, color = "gray") + 
        geom_point(aes(x = time_post_exp, y = value, color = data, shape = data), size = 3, alpha = 0.5) +
        scale_color_manual(values = c(color_study_A, color_study_B)) + 
        facet_grid(cols = vars(biomarker), rows = vars(antigen), scales = "free") + 
        theme_bw() + 
        theme( strip.background = element_blank(), legend.position = "top") + 
        labs(x = "Days after second dose", y = "Conc. of biomarker (%)", color = "Dataset", shape = "Dataset") + 
            david_theme_A() + 
        ggtitle("Biomarker data trajectories after second dose")

    p3 <- data_compare %>% filter(variable == "IC50hs") %>%
        ggplot() +
        geom_smooth(aes(x = time_post_exp, y = value), se = FALSE, color = "gray") + 
        geom_point(aes(x = time_post_exp, y = value, color = data, shape = data), size = 3, alpha = 0.5) +
        scale_color_manual(values = c(color_study_A, color_study_B)) + 
        facet_grid(cols = vars(biomarker), scales = "free") + 
        theme_bw() + 
        theme( strip.background = element_blank(), legend.position = "top") + 
        labs(x = "Days after second dose", y = "Log sVNT titre", color = "Dataset", shape = "Dataset") + 
            david_theme_A()

    require(patchwork)
    p0 <- p2 + p3 + plot_layout(widths = c(2, 1), guides = "collect") 
    p1 / p0  + plot_annotation(tag_levels = "A") & theme(legend.position='top') 
    ggsave(here::here("outputs", "figs", "data", "data_plot.png"), width = 10, height = 10)
}

#' A function to make table SM1 in the Appendix 
make_table1 <- function() {
    df_nih_meta <- readRDS(here::here("outputs", "data_clean", "nih_wu_s", "df_meta.RDS")) %>% mutate(data = "Study A (Calibration)")
    df_sid_meta <- readRDS(here::here("outputs", "data_clean", "sid_wu_s", "df_meta.RDS")) %>% mutate(data = "Study B (Validation)")

    library(table1) 
    library(flextable) 
    df_meta <- bind_rows(df_nih_meta, df_sid_meta)

    label(df_meta$vax_type)      <- "Vaccine type"
    label(df_meta$time_d1)      <- "Time since first dose (days)"
    label(df_meta$age_cat)    <- "Age (yrs)"



    tbl1 <- table1(~ vax_type + time_d1 + age_cat | data, data=df_meta)

    t1flex(tbl1) %>% 
    save_as_docx(path = here::here("outputs", "tables", "table1.docx"))
}