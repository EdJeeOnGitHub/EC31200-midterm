# This file replicates the main findings of 
# Dippel 2014
library(tidyverse) # data manipulation
library(haven) # reading stata files
library(testthat) # testing replication
library(fixest) # fixed effects/iv regression
library(broom) # tidy model extraction

############################################
## Data ####################################
############################################
# loading data
native_df <- read_dta("data/input/forcedcoexistence_webfinal.dta")
# subset of the data for 2000, used for most 
# tables and regressions
native_current_df <- native_df %>%
    filter(year == 2000) %>%
    rename(instrument_silver = intrument_silver) # fixing a typo in the orig


############################################
## Summary Stats ###########################
############################################
## Summary Stats
table_1 <- native_current_df %>%
    group_by(FC, HC) %>%
    summarise(
        mean_inc = mean(logpcinc),
        sd = sd(logpcinc),
        n = n()
    ) %>%
    pivot_wider(
        id_cols = c(
            FC,
            HC
        ), 
        names_from = c(
            FC
        ), 
        values_from = c(
            mean_inc:n
        ),
        names_glue = "FC{FC}_{.value}"
    )

# Some tests to ensure we're on the right track:
test_that("Table 1 Replicates", {
    # Forced Coexistence and Hist Cent {0,0}
    expect_equal(
        table_1[[1, "FC0_mean_inc"]],
        9.34,
        tolerance = 0.001)
    # Forced Coex and Hist Cent {0,1}
    expect_equal(
        table_1[[2, "FC0_mean_inc"]],
        9.58,
        tolerance = 0.001
    )
})


############################################
## OLS #####################################
############################################
# Table 2
### Create lists of controls and base OLS formula
base_ols_formula <- "logpcinc ~ FC + HC"
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

#'  Function to create OLS formulae
#'  
#' @param base_formula y and main independent variables of interest.
#' @param controls covariates we wish to control for, must be a string
create_ols_formula <- function(base_formula = base_ols_formula,
                               controls) {
    model_formula <- paste0(base_formula, " + ", paste(controls, collapse = " + "))
    return(as.formula(model_formula))
}

# Create formulae
ols_formulae <- list(
    "1", # base formula + intercept
    reservation_controls, 
    c(reservation_controls, tribe_controls),
    c(reservation_controls, tribe_controls, additional_reserve_controls),
    c(reservation_controls, tribe_controls, additional_reserve_controls, "0 | statenumber") # + state FE
) %>%
    map(~create_ols_formula(base_ols_formula,
                            .x))

# Fit OLS models
ols_models <- map(ols_formulae,
    ~feols(
        fml = .x,
        data = native_current_df,
        cluster = ~ eaid + statenumber
    ))


# Checking we have replicated the tables correctly
test_that("Table 2 OLS replications", {
    # read directly off table 2A
    dippel_FC_ols_coefs <- c(
        -0.358,
        -0.334,
        -0.364,
        -0.302,
        -0.291
    )
    # Estimated using our models
    # following pipe just extracts
    # coefs
    rep_FC_ols_coefs <- c(
        ols_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "FC") %>%
    map(select, estimate) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_FC_ols_coefs,
        dippel_FC_ols_coefs,
        expect_equal,
        tolerance = 0.001
    )
    # now t stats

    dippel_FC_ols_tstats <- c(
        -3.662,
        -4.090,
        -7.192,
        -4.913,
        -5.072
    )

    rep_FC_tstats <- c(
        ols_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "FC") %>%
    map(select, statistic) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_FC_tstats,
        dippel_FC_ols_tstats,
        expect_equal,
        tolerance = 0.001
    )

})

# creating dictionaries to clean up terms:
group_replacement <- list(
        "Reservation Controls" = reservation_controls,
        "Tribe Controls" = tribe_controls,
        "Extra Reservation Controls" = additional_reserve_controls,
        "Extra IV Controls" = additional_iv_controls
)
dict_replacment <- c(
        eaid = "Tribe",
        statenumber = "State",
        logpcinc = "log(per capita income)",
        FC = "Forced coexistence",
        HC = "Historical centralisation",
        instrument_silver = "Historical silver-mining",
        instrument_gold = "Historical gold-mining"
)
# aggregating models and saving to tex file
etable(
    ols_models,
    title = "OLS - Dippel Table 2, Panel A",
    drop = c("(Intercept)"),
    group = group_replacement,
    dict = dict_replacment,
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-2-ols.tex",
    digits = 3,
    coefstat = "tstat"
)

############################################
## FE ######################################
############################################
# Now panel B
# Create FE formulae
# We:
#   - Add "|eatribe" FE to the OLS formula,
#   - convert formulae to strings so we can remove brackets
#       and fix the case that already has state FE as the  syntax
#       is scuffed otherwise.
#   - convert back into formulae objects
fe_formulae <- ols_formulae %>%
    map(~update(.x, ~ . | eaid)) %>%
    map(deparse) %>%
    map(~str_remove_all(.x, "\\(|\\)") %>%paste(collapse = "")) %>%
    map(~str_replace_all(.x, "(?<=statenumber )\\|(?= eaid)", "+")) %>%
    map(as.formula)
# Fit FE models
fe_models <- map(fe_formulae,
    ~feols(
        fml = .x,
        data = native_current_df,
        cluster = ~eaid + statenumber
    ))



# Testing replication
test_that("Table 2 Panel B Replicates", {
    dippel_FC_fe_coefs <- c(
        -0.401,
        -0.318,
        -0.318,
        -0.274,
        -0.276
    )
    # Estimated using our models
    # following pipe just extracts
    # coefs
    rep_FC_fe_coefs <- c(
        fe_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "FC") %>%
    map(select, estimate) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_FC_fe_coefs,
        dippel_FC_fe_coefs,
        expect_equal,
        tolerance = 0.001
    )

})


etable(
    fe_models,
    title = "FE - Dippel Table 2, Panel B",
    group = group_replacement,
    dict = dict_replacment,
    coefstat = "tstat",
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-2-fe.tex",
    digits = 3,
    fitstat = ""
)

############################################
## FS/RF ###################################
############################################
##### First Stage Check
# Creating list of panel A first stage formulae
first_stage_A_formulae <- map(
    list(
        "1",
        reservation_controls,
        c(reservation_controls, tribe_controls),
        c(reservation_controls, tribe_controls, additional_reserve_controls),
        c(reservation_controls, tribe_controls, additional_reserve_controls, "0 |statenumber"),
        c(reservation_controls, tribe_controls, additional_reserve_controls, additional_iv_controls, "0 |statenumber")
    ),
    ~create_ols_formula(
        base_formula = "FC ~ instrument_gold + instrument_silver + HC",
        controls = .x
    )
)
# Same for panel B
first_stage_B_formulae <- map(
    list(
        "1",
        reservation_controls,
        c(reservation_controls, tribe_controls),
        c(reservation_controls, tribe_controls, additional_reserve_controls),
        c(reservation_controls, tribe_controls, additional_reserve_controls, "0 |statenumber"),
        c(reservation_controls, tribe_controls, additional_reserve_controls, additional_iv_controls, "0 |statenumber")
    ),
    ~create_ols_formula(
        base_formula = "logpcinc ~ instrument_gold + instrument_silver + HC",
        controls = .x
    )
)

# Fitting models
first_stage_A_models <- first_stage_A_formulae %>%
    map(~feols(
        fml = .x,
        data = native_current_df,
        cluster=~eaid + statenumber
    ))
first_stage_B_models <- first_stage_B_formulae %>%
    map(~feols(
        fml = .x,
        data = native_current_df,
        cluster=~eaid+statenumber
    ))
# Testing replication
test_that("Table IV  replicates ", {
    # read directly off table 4A
    dippel_fs_coefs <- c(
        0.018,
        0.017,
        0.015,
        0.018,
        0.023,
        0.030
    )
    # Estimated using our models
    # following pipe just extracts
    # coefs
    rep_fs_coefs <- c(
        first_stage_A_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "instrument_gold") %>%
    map(select, estimate) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_fs_coefs,
        dippel_fs_coefs,
        expect_equal,
        tolerance = 0.001
    )

})

etable(
    first_stage_A_models,
    title = "First Stage - Dippel Table 4, Panel A",
    group = group_replacement,
    drop = c("(Intercept)"),
    dict = dict_replacment,
    coefstat = "tstat",
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-4-A.tex",
    digits = 3
)
etable(
    first_stage_B_models,
    title = "First Stage - Dippel Table 4, Panel B",
    group = group_replacement,
    drop = c("(Intercept)"),
    dict = dict_replacment,
    coefstat = "tstat",
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-4-B.tex",
    digits = 3
)


############################################
## IV  #####################################
############################################

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


# Testing IV replication
test_that("Table V IV replications - two instruments", {
    # read directly off table 2A
    dippel_FC_iv_coefs <- c(
        -0.329,
        -0.304,
        -0.360,
        -0.316,
        -0.302,
        -0.403
    )
    # Estimated using our models
    # following pipe just extracts
    # coefs
    rep_FC_iv_coefs <- c(
        iv_2_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "fit_FC") %>%
    map(select, estimate) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_FC_iv_coefs,
        dippel_FC_iv_coefs,
        expect_equal,
        tolerance = 0.001
    )

})



etable(
    iv_2_models,
    title = "IV - Dippel Table 5, Panel A",
    group = group_replacement,
    drop = c("(Intercept)",
             "Historical centralisation"),
    dict = dict_replacment,
    coefstat = "tstat",
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-5-iv-A.tex",
    digits = 3,
    digits.stats = 3,
    fitstat = ~  ivf1 + ivwald1
)

## Now we fit panel B for the one instrument case.
## Aggregating into one variable
native_current_df <- native_current_df %>%
    mutate(instrument_precious = instrument_gold + instrument_silver)
# Mapping string formulae into actual formulae with precious metal instrument
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



test_that("Table V IV replications - one instruments", {
    # read directly off table 2A
    dippel_FC_iv_coefs <- c(
        -0.406,
        -0.371,
        -0.397,
        -0.350,
        -0.339,
        -0.443
    )
    # Estimated using our models
    # following pipe just extracts
    # coefs
    rep_FC_iv_coefs <- c(
        iv_2_models
    ) %>%
    map(tidy) %>%
    map(filter, term == "fit_FC") %>%
    map(select, estimate) %>%
    map(pull) %>%
    unlist()
    # compare each element
    map2(
        rep_FC_iv_coefs,
        dippel_FC_iv_coefs,
        expect_equal,
        tolerance = 0.001
    )

})

etable(
    iv_1_models,
    title = "IV - Dippel Table 5, Panel B",
    group = group_replacement,
    drop = c("(Intercept)",
             "Historical centralisation"),
    dict = dict_replacment,
    coefstat = "tstat",
    tex = TRUE,
    replace = TRUE,
    file = "data/output/table-5-iv-B.tex",
    digits = 3,
    digits.stats = 3,
    fitstat = ~  ivf1 + ivwald1
)
