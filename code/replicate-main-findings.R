library(tidyverse)
library(haven)
library(testthat)
native_df <- read_dta("data/input/forcedcoexistence_webfinal.dta")


## Summary Stats
table_1 <- native_df %>%
    group_by(FC, HC) %>%
    filter(year == 2000) %>%
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

native_current_df <- native_df %>%
    filter(year == 2000)
library(fixest)


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


ols_formulae_2_4 <- list(
    reservation_controls,
    c(reservation_controls, tribe_controls),
    c(reservation_controls, tribe_controls, additional_reserve_controls)
) %>%
    map(~create_ols_formula(base_ols_formula,
                            .x))


ols_models_1 <- feols(
    fml = as.formula(base_ols_formula),
    data = native_current_df,
    cluster = ~eatribe + statenumber
)

ols_models_2_4 <- map(ols_formulae_2_4,
    ~feols(
        fml = .x,
        data = native_current_df,
        cluster = ~ eatribe + statenumber
    ))

ols_models_2_4
ols_models_5 <- feols(
    fml = as.formula(paste0(paste0(deparse(ols_formulae_2_4[[3]]), collapse = ""), "| statenumber", collapse = "")),
    data = native_current_df,
    cluster = ~eatribe + statenumber
)
ols_models_5
native_current_df %>% 
    select(state, statenumber) %>%
    group_by(state) %>%
    distinct() %>%
    arrange(state)
    unique()

ols_models_2_4
ols_models_5

feols(ols_formula_3,
      native_current_df,
      cluster = ~eatribe+state)
model_fit <- feols(
    fml = ols_formula_2,
    data = native_current_df,
    cluster =  ~ eatribe + state 
)
model_fit
