# This file replicates the main findings of 
# Dippel 2014
library(tidyverse) # data manipulation
library(haven) # reading stata files
library(testthat) # testing replication
library(fixest) # fixed effects/iv regression
library(broom) # tidy model extraction
# loading data
native_df <- read_dta("data/input/forcedcoexistence_webfinal.dta")
# subset of the data for 2000, used for most 
# tables and regressions
native_current_df <- native_df %>%
    filter(year == 2000)


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

# Table 2




# Table 3


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

})

# creating dictionaries to clean up terms:
group_replacement <- list(
        "Reservation Controls" = reservation_controls,
        "Tribe Controls" = tribe_controls,
        "Extra Reservation Controls" = additional_reserve_controls
)
dict_replacment <- c(
        eaid = "Tribe",
        statenumber = "State",
        logpcinc = "log(per capita income)",
        FC = "Forced coexistence",
        HC = "Historical centralisation"
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
    digits = 3
)



