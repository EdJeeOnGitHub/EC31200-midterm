# This file extends the analysis of 
# Dippel 2014
library(tidyverse) # data manipulation
library(haven) # reading stata files
library(fixest) # fixed effects/iv regression
library(broom) # tidy model extraction
library(ivpack)
library(ggstance)

############################################
## Data ####################################
############################################
# loading data
native_df <- read_dta("data/input/forcedcoexistence_webfinal.dta")
# subset of the data for 2000, used for most 
# tables and regressions
native_current_df <- native_df %>%
    filter(year == 2000) %>%
    rename(instrument_silver = intrument_silver) %>% # fixing a typo in the orig
    mutate(instrument_precious = instrument_gold + instrument_silver)
reservation_controls <- c(
    "logpcinc_co",
    "logunempl_co",
    "logdist",
    "logruggedness",
    "logresarea_sqkm"
)
tribe_controls <- c(
    "ea_v5",
    "ea_v30",
    "ea_v32",
    "ea_v66"
)
additional_reserve_controls <- c(
    "logpop",
    "logpopsq",
    "popadultshare",
    "casino"
)

additional_iv_controls <- c(
    "homelandruggedness",
    "removal",
    "wgold_enviro",
    "wsilver_enviro"
)

# IV time
#' Create formulae for IV estimation
#'
#' @param controls covariates, strings.
#' @param fixed_effects which FE to control for, string. Leave NULL if unused.
#' @param instruments String of instruments to use
#'
#' @return IV formula
create_iv_formula <- function(controls,
                              fixed_effects = NULL,
                              instruments) {
    base_formula <- "logpcinc ~ HC"
    base_control <- paste0(base_formula, " + ", paste(controls, collapse = " + "))
    
    if (!is.null(fixed_effects)) {
        base_control_fe <- paste0(base_control, " | ", fixed_effects)
    } else {
        base_control_fe <- base_control
    }

    base_control_fe_iv <- paste0(base_control_fe, "| FC ~", instruments) 

    final_formula <- as.formula(base_control_fe_iv)
    return(final_formula)
}

# List of model formulae
iv_formulae_list <- list(
    "HC", # base formula + intercept
    c("HC", reservation_controls), 
    c("HC", reservation_controls, tribe_controls),
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls),
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls, "0 | statenumber"), # + state FE
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls, additional_iv_controls, "0 | statenumber")
) 

# mapping string list objects into actual formulae objects in R
iv_2_formulae <- iv_formulae_list %>%
    map(~create_iv_formula(
        controls = .x,
        instruments = "instrument_gold + instrument_silver"
    ))




# Fitting models
iv_2_models <- iv_2_formulae %>%
    map(~feols(
        fml = .x,
        data = native_current_df,
        cluster = ~eaid + statenumber
    ))  


iv_1_formulae <- iv_formulae_list %>%
    map(~create_iv_formula(
        controls = .x,
        instruments = "instrument_precious"
    ))


iv_1_models <- iv_1_formulae %>%
    map(~feols(
        fml = .x,
        data = native_current_df,
        cluster = ~eaid + statenumber
    ))

create_iv_reg_formula <- function(dependent_var,
                                  endog_var,
                                  controls,
                                  instruments){
    base_formula <- paste0(dependent_var, " ~ ", endog_var)
    iv_formula <- paste0(
        base_formula,
        " + ",
        paste(controls, collapse = " + "),
        "| ", 
        paste(instruments, collapse = " + "), " + ",
        paste(controls, collapse = " + "))
    final_formula <- as.formula(iv_formula)
    return(final_formula)
}

ivreg_formulae_list <- list(
    "HC", # base formula + intercept
    c("HC", reservation_controls), 
    c("HC", reservation_controls, tribe_controls),
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls),
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls, "factor(statenumber)"), # + state FE
    c("HC", reservation_controls, tribe_controls, additional_reserve_controls, additional_iv_controls, "factor(statenumber)")
) 


ivreg_one_instrument_formulae <- ivreg_formulae_list %>%
    map(~create_iv_reg_formula(
        dependent_var = "logpcinc",
        endog_var = "FC",
        controls = .x,
        instruments = c("instrument_precious")
    ))

ivreg_one_instrument_models <- ivreg_one_instrument_formulae %>%
    map(~ivreg(
        formula = .x,
        data = native_current_df,
        x = TRUE
    ))




one_instrument_conf_interval <- ivreg_one_instrument_models %>%
    map_dfr(anderson.rubin.ci) %>%
    mutate(ar.conf.low = str_extract(confidence.interval, "(?<=\\[).*(?=,)")%>%as.numeric,
           ar.conf.high = str_extract(confidence.interval, "(?<=,).*(?=\\])")%>%as.numeric) %>%
    mutate(model = 1:n()) %>%
    select(-confidence.interval)

iv_1_confidence_df <- inner_join(
    
    iv_1_models %>%
        imap_dfr(~tidy(.x, conf.int =  TRUE) %>% mutate(model = .y)) %>%
        filter(term == "fit_FC") %>%
        rename(trad.conf.high = conf.high,
               trad.conf.low = conf.low),
    one_instrument_conf_interval,
    by = "model"
)
iv_1_long_confidence_df <- iv_1_confidence_df %>%
    pivot_longer(cols = contains("conf"),
    names_to = c("type", "hilo"),
    names_pattern = "(.*)\\.conf\\.(.*)") %>%
    spread(hilo, value)

anderson_rubin_plot <- iv_1_long_confidence_df %>%
    mutate(model = factor(model, levels = 6:1)) %>%
    ggplot(aes(x = estimate,
               xmin = low,
               xmax = high,
               y = model,
               colour = type)) +
    geom_pointrangeh(position = position_dodge2v(0.5)) +
    geom_vline(xintercept = 0,
               linetype = "longdash") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_discrete(labels = c("Anderson-Rubin", "Conventional")) +
    labs(y = "Model",x = "Estimate")
ggsave(
    plot = anderson_rubin_plot,
    filename = "data/output/anderson-rubin-plot.png",
    width = 8,
    height = 6)