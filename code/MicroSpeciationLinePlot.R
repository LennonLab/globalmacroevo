## Plotting Microbial Speciation against Extinction Rates
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

# import Rohde & Muller extinction events data
# Originally from MicroSpeciation-Logistic3

# Say x is diversity and y is extinction intensity for an era lasting 3 years
# Current dataset has data like the following for an era lasting 3 My:
# Time  Diversity   ExtinctionIntensity
# t     x           y 
# t+1   x           y
# t+2   x           y

# I'm altering it to look like this to spread out extinction over
# the entire length of the era:

# Time  ExtinctionIntensity
# t     y/3 
# t+1   y/3
# t+1   y/3

ext_file <- paste(data_dir, "Extinction_Rohde_Muller_raw.csv", sep = "")
ext <- read.csv(ext_file)

avg_ext <- c()

# Average out exitinction percentage over entire time period
for (i in 1:nrow(ext)){
  
  if (i == 1){
    # First entry in the dataframe starts a paleontologic era
    era <- c(ext$ExtinctionIntensity[i])
    
    # For all entries after the first: 
  } else {
    # I assume that when either ExtinctionIntensity or Diversity changes,
    # then a new era has begun, otherwise, an era continues
    if (ext$Diversity[i] == ext$Diversity[i-1] & ext$ExtinctionIntensity[i] == ext$ExtinctionIntensity[i-1]) {
      era <- c(era, ext$ExtinctionIntensity[i])
      if (i == nrow(ext)) { # if the last entry is its own era
        
        # integer(x) makes a vector of zeros of length x
        # the mean of extinction intensity is evenly distributed over an era's duration
        # avg_ext <- c(avg_ext, integer(length(era)) + mean(era)/length(era))
        a <- integer(length(era))
        a[length(a)] <- mean(era)
        avg_ext <- c(avg_ext, a)
      }
      
    } else { # if a new era is beginning
      # Add previous era into vector
      a <- integer(length(era))
      a[length(a)] <- mean(era)
      avg_ext <- c(avg_ext, a)
      era <- c(ext$ExtinctionIntensity[i]) # begin recording new era
      
      if (i == nrow(ext)) { # if the last entry is part of another era
        a <- integer(length(era))
        a[length(a)] <- mean(era)
        avg_ext <- c(avg_ext, a)
      }
    }
  }
}

# Add back into dataframe
ext <- data.frame(ext, AverageExtinction = avg_ext) 
ext$time <- as.character(4000 - ext$Mya)

names(avg_ext) <- ext$time
avg_ext <- avg_ext[order(names(avg_ext))]/100 

## simulations

# timestep function - discrete time speciation
timestep <- function(S, r, ep = 0){
  return(S * exp(r)* (1 - ep))
}

lambda <- c(0.020, 0.025, 0.031) # speciation
epsilon <- 0.7 # some extinction rates
S1 <- c(1) # initial species levels 4Ga - first parameter set
S2 <- c(1) # initial species levels 4Ga - second parameter set
S3 <- c(1) # initial species levels 4Ga - third parameter set
EeStart <- as.integer(names(avg_ext[1])) # first year with extinction data

# Simulate, but replace random background extinction levels 
# with "known" extinction rates for macroorganisms for "recent" years
for(t in 2:4000) {
  if (t < EeStart){ # before data, make extinction be random
    ep <- rnorm(n = 1, mean = epsilon, sd = 0.1) # give me an epsilon near the actual epsilon for all simulations per time step
    ep <- ifelse(ep<0,0.01,ep) # if epsilon is ever less than 0, make it 0.1
    ep <- ifelse(ep>1,0.9,ep) # don't want epsilon above 0.9
    mu <- lambda * ep 
    r <- lambda - mu # diversification rate
    
    S1[t] <- timestep(S = S1[t-1], r = r[1]) # discrete exponential growth, pure birth
    S2[t] <- timestep(S = S2[t-1], r = r[2])
    S3[t] <- timestep(S = S3[t-1], r = r[3]) # discrete exponential growth, pure birth
    
    
  } else { # when extinction data starts, use it in place of random data
    r <- lambda # diversification rate is just lambda
    # ep <- avg_ext[t - EeStart + 1]
    if (t < 4000){
      ep <- avg_ext[t - EeStart + 1]
    } else {
      ep <- median(avg_ext)
    }
    
    S1[t] <- timestep(S = S1[t-1], r = r[1], ep = ep) # ep % of population dies
    S2[t] <- timestep(S = S2[t-1], r = r[2], ep = ep)
    S3[t] <- timestep(S = S3[t-1], r = r[3], ep = ep)
  }
}
df <- data.frame(logrichness = c(log10(S1), log10(S2), log10(S3)), time = c(1:4000, 1:4000, 1:4000), lambda = c(double(4000) + lambda[1], double(4000) + lambda[2], double(4000) + lambda[3]))
(p1 <- ggplot(data = df, aes(y = logrichness, x = time, group = as.character(lambda))) + 
  geom_line(aes(linetype = as.character(lambda))) +
  scale_x_continuous("Time (Myr)", expand = c(0,0)) +
  scale_y_continuous("Log(Richness)", breaks = seq(0,32,4), expand = c(0,0)) +
  scale_linetype_discrete(expression(lambda)) +
  geom_vline(xintercept = EeStart, linetype = "dotted", color = "red") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        axis.line = element_line(colour = "black")
))

fig_dir <- paste(figure_dir, "ExtinctionEvents_lambdalineplots.png", sep = "")
ggsave(plot = p1, filename = fig_dir, width = 7, height = 5)
                                                                                             