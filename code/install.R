# Installation script
install <- FALSE
package_vector <-
c(
    "tidyverse",
    "ivpack",
    "fixest",
    "broom",
    "ggstance",
    "maps"
    "tidygeocoder",
    "ggrepel",
    "spdep",
    "kableExtra"
)                             
if (install == TRUE){
    lapply(package_vector, install.packages)
    devtools::install_github("https://github.com/rgiordan/zaminfluence/", 
                               ref="master", 
                               subdir="zaminfluence", 
                               force=TRUE) 
}