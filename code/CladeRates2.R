# Simulate variable rates among clades
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

species <- list("0.1" = double(4000), "0.5" = double(4000), "0.9" = double(4000))
species$`0.1`[1] = 1.0
S_total <- double(4000)
S_total[1] <- 1.0
lambda <- 0.025
epsilon <- c(0.1, 0.5, 0.9)
mu <- epsilon * lambda
r <- lambda - mu
timestep <- function(species, S_total, i){
  species$`0.1`[i] <-  0.05 * species$`0.5`[i-1] * (r[2]) +
    0.97 * species$`0.1`[i-1] * (r[1]) +
    0.05 * species$`0.9`[i-1] * (r[3]) +
    0.95* species$`0.1`[i-1]
  species$`0.5`[i] <-  0.90 * species$`0.5`[i-1] * (r[2]) +
    0.0 * species$`0.1`[i-1] * (r[1]) +
    0.05 * species$`0.9`[i-1] * (r[3]) +
    0.95 *species$`0.5`[i-1] +
    0.04 * species$`0.1`[i-1]
  species$`0.9`[i] <-  0.05 * species$`0.5`[i-1] * (r[2]) +
    0 * species$`0.1`[i-1] * (r[1]) +
    0.90 * species$`0.9`[i-1] * (r[3]) +
    species$`0.9`[i-1] +
    0.05 *species$`0.5`[i-1] +
    0.01 * species$`0.1`[i-1]

  S_total[i] <- species$`0.1`[i] +
                 species$`0.5`[i] +
                 species$`0.9`[i]
  return(list(species = species, S_total = S_total))
}
# timestep <- function(species, S_total, i){
#   species$`0.1`[i] <-  0.015 * species$`0.5`[i-1] * (r[2]) + 
#     0.97 * species$`0.1`[i-1] * (r[1]) +
#     0.02 * species$`0.9`[i-1] * (r[3]) +
#    species$`0.1`[i-1]
#   species$`0.5`[i] <-  0.97 * species$`0.5`[i-1] * (r[2]) + 
#     0.0 * species$`0.1`[i-1] * (r[1]) +
#     0.04 * species$`0.9`[i-1] * (r[3]) +
#     species$`0.5`[i-1] 
#   species$`0.9`[i] <-  0.015 * species$`0.5`[i-1] * (r[2]) + 
#     0 * species$`0.1`[i-1] * (r[1]) +
#     0.94 * species$`0.9`[i-1] * (r[3]) + 
#     species$`0.9`[i-1]
#   
#   S_total[i] <- species$`0.1`[i] + 
#     species$`0.5`[i] +
#     species$`0.9`[i]
#   return(list(species = species, S_total = S_total))
# }


for (i in 2:4000){
  a <- timestep(species, S_total, i)
  species <- a$species
  S_total <- a$S_total
}

df <- data.frame(logrichness = c(log10(species$`0.1`), log10(species$`0.5`), log10(species$`0.9`), log10(S_total)), 
                 time = c(1:4000, 1:4000,1:4000, 1:4000),
                 label = c(as.character(double(4000) + 0.1), as.character(double(4000) + 0.5), as.character(double(4000) + 0.9), paste(character(4000), "Total", sep = "")))
ggplot(df, aes(y = logrichness, x = time, color = as.character(label))) +
  geom_line() +
  scale_color_viridis("adad", discrete = T)
