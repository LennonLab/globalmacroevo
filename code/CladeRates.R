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

epsilon <- seq(0.1,0.9, by = 0.05)
species <- matrix(double(4000 * length(epsilon)),nrow = 4000, ncol = length(epsilon))
colnames(species) <-  as.character(epsilon)
species[1, "0.7"] = 1
S_total <- double(4000)
S_total[1] <- 1.0
lambda <- 0.015
mu <- epsilon * lambda
r <- lambda - mu
timestep <- function(ep, clades, ind){
  ep <- as.character(ep)
  i <- match(ep, epsilon) 
  if (i == 1){
    ep1 <- as.character(epsilon[i+1])
    ep2 <- as.character(epsilon[i+2])
    clades[ind, ep] <- 0.97 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      clades[ind-1, ep]
  } else if (i == 2){
    ep1_ <- as.character(epsilon[i-1])
    ep1 <- as.character(epsilon[i+1])
    ep2 <- as.character(epsilon[i+2])
    clades[ind, ep] <- 0.95 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      clades[ind-1, ep]
  } else if (epsilon[i] == epsilon[length(epsilon) -1]){
    ep1_ <- as.character(epsilon[i-1])
    ep2_ <- as.character(epsilon[i-2])
    ep1 <- as.character(epsilon[i+1])
    clades[ind, ep] <- 0.95 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1, ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  } else if (epsilon[i] == epsilon[length(epsilon)]){	
    ep1_ <- as.character(epsilon[i-1])
    ep2_ <- as.character(epsilon[i-2])
    clades[ind, ep] <- 0.97 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1,ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  } else{
    ep1_ <- as.character(epsilon[i-1])
    ep2_ <- as.character(epsilon[i-2])
    ep1 <- as.character(epsilon[i+1])
    ep2 <- as.character(epsilon[i+2])
    clades[ind, ep] <- 0.94 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1,ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  }
  return(clades)
}


for (i in 2:4000){
  for (ep in epsilon){
    species <- timestep(ep = ep, clades = species, ind = i)
    S_total[i] <- S_total[i] + species[i, as.character(ep)]
  }
  
}

df <- as.data.frame(species)
df$time <- 1:4000
species_df <- gather(df, key = epsilon, value = Richness, `0.1`:`0.9`)
species_df$Richness <- ifelse(species_df$Richness == 0, 10^-20, species_df$Richness)
(p1<- ggplot(species_df, aes(x = time, y = Richness, color = as.double(epsilon))) + 
  geom_line(aes(group = as.double(epsilon)), size = 0.3) +
  scale_x_continuous("Time (Myr)", breaks = c(0, 1000, 2000, 3000, 4000)) +
  scale_y_log10(limits = c(1, S_total[4000]), breaks = 10^(seq(1,14, by = 2))) +
  stat_function(fun = function(time) log10(S_total[time]), color = "black", linetype = "dashed") +
  scale_color_viridis(expression(epsilon), option = "C") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
)

fig_dir <- paste(figure_dir, "VariableCladeRates.png", sep = "")
ggsave(plot = p1, filename = fig_dir, width = 7, height = 5)
