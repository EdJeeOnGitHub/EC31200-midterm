# This file extends the analysis of 
# Dippel 2014
library(tidyverse) # data manipulation
library(haven) # reading stata files
library(fixest) # fixed effects/iv regression
library(broom) # tidy model extraction
library(ivpack)
library(ggstance)
library(zaminfluence)
library(spdep)
library(sf)
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


###
binary_native_df <-  native_current_df %>% 
    select(statenumber, FC, instrument_precious, logpcinc) %>%
    group_by(statenumber) %>%
    mutate(mean_d = median(FC),
           mean_z = median(instrument_precious),
           d = as.numeric(FC > mean_d),
           z = as.numeric(instrument_precious > mean_z))

find_proportions <- function(df){
    count_df <- df %>%
    group_by(z, d) %>%
    summarise(n = n())
    
    p_df <- count_df %>%
        ungroup() %>%
        group_by(z) %>%
        mutate(pr = n / sum(n)) 

    p_always <- p_df %>%
        filter(z == 0 & d == 1) %>%
        summarise(p_always = sum(pr)) %>% 
        select(p_always) %>%
        pull()
    p_never <- p_df %>%
        filter(z == 1 & d == 0) %>%
        summarise(p_never = sum(pr)) %>%
        select(p_never) %>%
        pull()

    p_complier <- 1 - p_never - p_always
    n_complier <- nrow(df) * p_complier
    return(list(
        p_always = p_always,
        p_never = p_never,
        p_complier = p_complier,
        n_complier = n_complier
    ))    
}

p_c_all <- binary_native_df %>%
    find_proportions()
p_c_all
min_y <- min(binary_native_df$logpcinc)
max_y <- max(binary_native_df$logpcinc)

dens_df <- binary_native_df %>%
    group_by(z, d) %>%
    summarise(dens = list(
        density(
            logpcinc,
            from = min_y,
            to = max_y,
            n = 100
        )$y
    )) %>%
    pivot_wider(values_from = dens, names_from = c(z, d)) %>%
    mutate(g_c_0 = list(unlist(`0_0`) *(p_c_all$p_complier + p_c_all$p_never)/p_c_all$p_complier - unlist(`1_0`)*(p_c_all$p_never/p_c_all$p_complier)),
           g_c_1 = list(unlist(`1_1`) *(p_c_all$p_complier + p_c_all$p_always)/p_c_all$p_complier - unlist(`0_1`)*(p_c_all$p_always/p_c_all$p_complier)),
           g_n = `1_0`,
           g_a = `0_1`,
           y = list(seq(from = min_y, to = max_y, length.out = 100)))  %>%
    select(contains("g"), y) %>%
    gather(variable, value, -y) %>%
    unnest() %>%
    mutate(variable_type = case_when(
        variable == "g_a" ~ "Always/Never",
        variable == "g_n" ~ "Always/Never",
        variable == "g_c_0" ~ "Complier", 
        variable == "g_c_1" ~ "Complier"
    )) 
dens_df
 
dens_df %>%
    group_by(variable) %>%
    summarise(mean_val = sum(y*value)*(max_y - min_y)/600%>%round(2)) %>%
    ungroup() %>%
    mutate(variable = case_when(
        variable == "g_a" ~ "Always Takers",
        variable == "g_n" ~ "Never Takers",
        variable == "g_c_0" ~ "Complier Y_0", 
        variable == "g_c_1" ~ "Complier Y_1"
    )) 

density_plot <- dens_df %>%
        mutate(variable = case_when(
        variable == "g_a" ~ "Always Takers",
        variable == "g_n" ~ "Never Takers",
        variable == "g_c_0" ~ "Complier Y_0", 
        variable == "g_c_1" ~ "Complier Y_1"
    )) %>%  
    mutate(variable = factor(variable, levels = c("Always Takers", "Never Takers", "Complier Y_0", "Complier Y_1"))) %>%
    ggplot(aes(x = y, y = value, color = variable, fill = variable)) +
    geom_line(size = 1) +
    theme_bw() + 
    theme(legend.position =  "bottom") + 
    labs(title = "Counterfactual/Observed Densities") +
    facet_wrap(~variable_type)

ggsave(
    plot = density_plot,
    filename = "data/output/density-plot.png",
    width = 8,
    height = 6)


######### Run AMIP ########
iv_infl <- ComputeModelInfluence(ivreg_one_instrument_models[[6]], se_group = native_current_df$statenumber)
grad_df <- GetTargetRegressorGrads(iv_infl, "FC")
influence_dfs <- SortAndAccumulate(grad_df)

change_df <- GetRegressionTargetChange(influence_dfs, "num_removed") %>%
    na.omit() 
change_df

######### Running Spatial Noise #############
run_geocode <- FALSE
if (run_geocode == TRUE) {

    library(tidygeocoder)

    geo_native_df <- native_current_df %>%
        mutate(clean_geo_name = str_remove_all(geo_name, " and Off-Reservation Trust Land"),
            clean_geo_name = paste0(clean_geo_name, ", USA")) %>%
        geocode(clean_geo_name, method = 'osm', lat = latitude , long = longitude)
    write_csv(geo_native_df, "data/output/geo-native-df.csv")
} else {
    geo_native_df <- read_csv("data/output/geo-native-df.csv")
}




ggplot(geo_native_df, aes(longitude, latitude), color = "grey99") +
  borders("state") + geom_point() +
#   geom_label_repel(aes(label = clean_geo_name)) +
  theme_void()
#
geo_native_sf <- geo_native_df %>%
    filter(!is.na(latitude)) 
coordinates(geo_native_sf) <- ~ longitude + latitude

neighbours <- knearneigh(coordinates(geo_native_sf), 3, longlat = TRUE) %>%
    knn2nb()

neighbours_list <- nb2listw(neighbours)

moran_sim <- moran.mc(
    geo_native_sf$logpcinc,
    neighbours_list,
    nsim = 100000
)
plot(moran_sim)
moran_df <- tibble(x = moran_sim$res)
moran_plot <- moran_df %>%
    ggplot(aes( x = x)) +
    geom_histogram(bins = 60, colour = "white", fill = "hotpink") +
    geom_vline(xintercept = moran_sim$statistic,
               linetype = "longdash") +
    theme_bw() +
    labs(title = "Distribution of Permuted Test Statistics and Realised Draw",
         x = "Estimate")
moran_plot
ggsave(plot = moran_plot,
       filename = "data/output/moran-plot.png",
       width = 8,
       height = 6)



library(plgp)
clean_geo_native_df <- geo_native_df %>%
    filter(!is.na(longitude)) 
dist_matrix <- clean_geo_native_df %>%
    select(longitude, latitude) %>% 
    distance()
eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-D) + diag(eps, nrow(clean_geo_native_df))
sim_y_draws <- mvtnorm::rmvnorm(1000, sigma = Sigma) %>%
    t() %>%
    as_tibble()



fit_sim <- function(y_draws){

    sim_data <- bind_cols(
        clean_geo_native_df,
        y_draws
    )
    sim_fit <- feols(
        fml = value ~ HC  + logpcinc_co + logunempl_co + logdist + logruggedness + 
    logresarea_sqkm + ea_v5 + ea_v30 + ea_v32 + ea_v66 + logpop + 
    logpopsq + popadultshare + casino + homelandruggedness + 
    removal + wgold_enviro + wsilver_enviro + 0 | statenumber | 
    FC ~ instrument_precious,
    data = sim_data
    ) %>%
    tidy()
    return(sim_fit)
}

nested_sim_y_draws <- sim_y_draws %>%
    gather(variable, value) %>%
    group_by(variable) %>%
    nest()
simulation_fits <- nested_sim_y_draws %>%
    mutate(sim_fit = map(data, ~fit_sim(y_draws =  .x)))

simulation_fits %>%
    unnest(sim_fit) %>%
    filter(term == "fit_FC") %>%
    ggplot(aes(x = statistic)) + 
    geom_histogram() +
    geom_vline(xintercept = realised_t_stat)



realised_t_stat <- iv_1_models[[6]] %>%
    tidy() %>%
    filter(term == "fit_FC") %>%
    select(statistic) %>%
    pull()



realised_pval <- iv_1_models[[6]] %>%
    tidy() %>%
    filter(term == "fit_FC") %>%
    select(p.value) %>%
    pull()    
simulation_fits %>%
    unnest(sim_fit) %>%
    filter(term == "fit_FC") %>%
    ggplot(aes(x = p.value)) + 
    geom_histogram() 

qq_plot <- simulation_fits %>%
    unnest(sim_fit) %>%
    filter(term == "fit_FC") %>%
    ggplot(aes(sample = p.value)) + 
    geom_qq() +
    geom_qq_line() + 
    theme_bw() +
    labs(title = "QQ-plot of Dependent Variable Spatial Noise P-Values")
ggsave(plot = qq_plot,
       filename = "data/output/qq-plot.png",
       width = 8,
       height = 6) 
