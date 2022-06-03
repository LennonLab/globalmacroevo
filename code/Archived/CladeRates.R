## Simulate variable rates among clades
## Ford Fishman

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
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

# Set up parameters
time <- 4000
epsilon <- seq(0.1,0.9, by = 0.05) # possible range of relative extinction rates (epsilon)
# set up matrix with columns being clade extinction rate and rows representing species over timr
species <- matrix(double(time * length(epsilon)),nrow = time, ncol = length(epsilon)) 
colnames(species) <-  as.character(epsilon)
# starting species has the following relative extinction rate:
species[1, "0.7"] = 1
S_total <- double(time)
S_total[1] <- 1.0 # only 1 starting species
lambda <- 0.015 # speciation/my, set for the entire simulation
mu <- epsilon * lambda # extinction/my
r <- lambda - mu # diversification/my

# function similating diversification process
# Clades at a given level of relative extinction speciate into clades with slighly lower or higher extinction
# Up to 6% of a clades newly speciated members change their relative extinction rate
# Arguments:
# epsilon - the given relative extinction rate of the current clade
# clades - the richness levels of all clades organized by time and clade
# ind - the current time

timestep <- function(ep, clades, ind){
  ep <- as.character(ep) # what epsilon does this clade have
  i <- match(ep, epsilon) # index of that epsilon value in the vector
  # needed to know which other possible values of epsilon are nearby 
  
  # if this is the lowest value of epsilon
  if (i == 1){
    
    ep1 <- as.character(epsilon[i+1]) # next higher level of epsilon
    ep2 <- as.character(epsilon[i+2]) # 2 levels up for epsilon
    
    clades[ind, ep] <- 0.97 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      clades[ind-1, ep]
  # if this is the second lowest value of epsilon
  } else if (i == 2){
    
    ep1_ <- as.character(epsilon[i-1]) # next lower level of epsilon
    ep1 <- as.character(epsilon[i+1]) # next higher level of epsilon
    ep2 <- as.character(epsilon[i+2]) # 2 levels up for epsilon
    
    clades[ind, ep] <- 0.95 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      clades[ind-1, ep]
  # if this is the second highest level of epsilon
  } else if (epsilon[i] == epsilon[length(epsilon) -1]){
    
    ep1_ <- as.character(epsilon[i-1]) # next lower level of epsilon
    ep2_ <- as.character(epsilon[i-2]) # 2 levels lower for epsilon
    ep1 <- as.character(epsilon[i+1]) # next higher level of epsilon
    
    clades[ind, ep] <- 0.95 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1, ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  # if this is the highest level of epsilon  
  } else if (epsilon[i] == epsilon[length(epsilon)]){	
    
    ep1_ <- as.character(epsilon[i-1]) # next lower level of epsilon
    ep2_ <- as.character(epsilon[i-2]) # 2 levels lower for epsilon
    
    clades[ind, ep] <- 0.97 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1,ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  # if this is any other value of epsilon
  } else{
    ep1_ <- as.character(epsilon[i-1]) # next lower level of epsilon
    ep2_ <- as.character(epsilon[i-2]) # 2 levels lower for epsilon
    ep1 <- as.character(epsilon[i+1]) # next higher level of epsilon
    ep2 <- as.character(epsilon[i+2]) # 2 levels up for epsilon
    
    clades[ind, ep] <- 0.94 * clades[ind-1, ep] * (r[i]) + 
      0.02 * clades[ind-1, ep1] * (r[i+1]) +
      0.01 * clades[ind-1, ep2] * (r[i+2]) +
      0.02 * clades[ind-1, ep1_] * (r[i-1]) +
      0.01 * clades[ind-1,ep2_] * (r[i-2]) +
      clades[ind-1, ep]
  }
  return(clades)
}

# run the function 
for (i in 2:time){
  for (ep in epsilon){
    species <- timestep(ep = ep, clades = species, ind = i)
  }
  
}

S_total <- rev(apply(species,MARGIN = 1,FUN = sum)) # total up the species across clades, put into reverse chronological order
df <- as.data.frame(species)
df$time <- time:1 # want years to start at 4000 mya
species_df <- gather(df, key = epsilon, value = Richness, `0.1`:`0.9`)
species_df$Richness <- ifelse(species_df$Richness == 0, 10^-20, species_df$Richness)

(p1<- ggplot(species_df, aes(x = time, y = Richness, color = as.double(epsilon))) + 
  geom_line(aes(group = as.double(epsilon)), size = 0.5) +
  scale_x_reverse("Time (Mya)", breaks = c(0, 1000, 2000, 3000, 4000)) +
  scale_y_log10(limits = c(1, S_total[1]), breaks = 10^(seq(1,14, by = 2))) +
  stat_function(fun = function(time) log10(S_total[time]), color = "black", linetype = "dashed") + # total number of species
  scale_color_viridis(expression(epsilon), option = "C") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 
)

fig_dir <- paste(figure_dir, "VariableCladeRates.png", sep = "")
ggsave(plot = p1, filename = fig_dir, width = 7, height = 5)
