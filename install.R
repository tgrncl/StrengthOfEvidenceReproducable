# Installing the required R packages

list.of.packages <- c("plyr", "dplyr",            # plyr and dplyr are used for data wrangling
                      "pwr", "esc", "compute.es", # pwr, esc, and comput.es compute effect size and power transformations
                      "car",                      # car is used for utility functions, such as outlier identification
                      "ggplot2",                  # ggplot2 is the main graphic engine, 
                      "ggExtra", "ggridges",      # ggExtra for marginal plots and ggridges for power ridge plots
                      "viridis",                  # viridis makes the standard color palette
                      "scico",                    # scico offers further palettes for a bipartitioned heatmap
                      "weights",                  # used for rounding w/o leading zero
                      "shiny")                    # shiny is responsible for the Web hosting
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
