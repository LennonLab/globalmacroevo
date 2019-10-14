## Simulate variable rates among clades
## Ford Fishman

## PURPOSE

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/Documents/LennonLab/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)

# Load packages
require("png")
require("grid")
require("tidyr")
require("ggplot2")
require("viridis")

species <- list("0.1" = c(1), "0.5" = c(0), "0.9" = c(0))
S_total <- c(1)
lambda <- 0.015
epsilon <- c(0.1, 0.5, 0.9)
mu <- epsilon * lambda
r <- lambda - mu
timestep <- function(species, S_total){
  s0.5 <-  0.30 * species$`0.5`[length(species$`0.5`)] * (r[2]) + 
    0.30 * species$`0.1`[length(species$`0.1`)] * (r[1]) +
    0.008 * species$`0.9`[length(species$`0.9`)] * (r[3]) + species$`0.5`[length(species$`0.5`)]
  s0.1 <-  0.30 * species$`0.5`[length(species$`0.5`)] * (r[2]) + 
    0.20 * species$`0.1`[length(species$`0.1`)] * (r[1]) +
    0.002 * species$`0.9`[length(species$`0.9`)] * (r[3]) + species$`0.1`[length(species$`0.1`)]
  s0.9 <-  0.40 * species$`0.5`[length(species$`0.5`)] * (r[2]) + 
    0.50 * species$`0.1`[length(species$`0.1`)] * (r[1]) +
    0.99 * species$`0.9`[length(species$`0.9`)] * (r[3]) + species$`0.9`[length(species$`0.9`)]
  species$`0.1` <- c(species$`0.1`, s0.1)
  species$`0.5` <- c(species$`0.5`, s0.5)
  species$`0.9` <- c(species$`0.9`, s0.9)
  S_total <- c(S_total,
               species$`0.1`[length(species$`0.1`)] + 
                 species$`0.5`[length(species$`0.5`)] +
                 species$`0.9`[length(species$`0.9`)])
  return(list(species = species, S_total = S_total))
}

  

for (i in 2:4000){
  a <- timestep(species, S_total)
  species <- a$species
  S_total <- a$S_total
}
S_total[4000]
head(species$`0.1`)
ggplot(NULL, aes(x = 1:4000)) +
  stat_function(fun = function(time) species$`0.9`[time], color = "black", linetype = "dashed") +
  stat_function(fun = function(time) species$`0.5`[time], color = "red", linetype = "dashed") +
  stat_function(fun = function(time) species$`0.1`[time], color = "blue", linetype = "dashed") +
  scale_y_log10()

  