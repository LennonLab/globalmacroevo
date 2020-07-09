## Modeling impacts of mass extiction on microbial diversity
## Does not use Sepkoski intensities
## Ford Fishman
library(scales)
library(MASS)
library(grid)
library(ggplot2)
library(tidyr)
library(png)

# Load packages and functions, set environment
source("~/GitHub/MicroSpeciation/code/functions.R")

# Mass extinction events
GOE <- 2450
OrdovicianSilurian <- 445
LateDevonian <- 375
PermianTriassic <- 252
TriassicJurassic <- 201
CretaceousPaleogene <- 66

# timestep function - discrete time speciation
timestep <- function(S, r){
  return(S * exp(r))
}

EE <- c(GOE, OrdovicianSilurian, LateDevonian, PermianTriassic,TriassicJurassic, CretaceousPaleogene)

lambda <- 0.015
epsilon <- rep(0.5,4)
S1 <- c(1) # initial species levels 4Ga - first parameter set
S2 <- c(1) # initial species levels 4Ga - second parameter set
S3 <- c(1) # initial species levels 4Ga - third parameter set
S4 <- c(1) # initial species levels 4Ga - fourth parameter set

for(t in 2:4000) {
  mya <- 4001 - t # how many millions years ago is it?
  
  ep <- rnorm(n = length(epsilon), mean = epsilon, sd = 0.3) # give an epsilon near the actual epsilon for all simulations per time step
  ep <- ifelse(ep<0,0,ep) # if epsilon is ever less than 0, make it 0
  # ep <- ifelse(ep>1,0.9,ep) # don't want epsilon above 0.
  
  if (mya %in% EE) {
    
    
    
    p <- c(0, 0.1, 0.5, 0.9) # proportion of taxa that go extinct
    # ep <- 0 # turn off background extinction
    
  } else {

    p <- rep(0,4) 
  }
  
  mu <- lambda * ep 
  r <- lambda - mu # diversification rate
  S1[t] <- S1[t-1]*exp(r[1])*(1-p[1])
  S2[t] <- S2[t-1]*exp(r[2])*(1-p[2])
  S3[t] <- S3[t-1]*exp(r[3])*(1-p[3])
  S4[t] <- S4[t-1]*exp(r[4])*(1-p[4])
}

df1 <- data.frame(zero = S1, ten = S2, fifty = S3, ninety = S4, mya=4000:1)

df2 <- gather(df1, key = "intensity", value = "richness",-mya)

# library("scales")

ggplot(df2, aes(x = mya, y = richness, group = intensity)) +  
  geom_line(aes(linetype = intensity)) +
  scale_y_log10("Taxon Diversity", expand = c(0,0)) +
  scale_x_continuous("Millions of Years Ago", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 10)) +
  geom_vline(xintercept = EE, size = 4, alpha = 0.4) +
  scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed","dotdash"),
                        limits = c("zero", "ten", "fifty","ninety"), 
                        labels = c("0%", "10%","50%", "90%")) +
  annotation_logticks(sides = "lb") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        # axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt")
  )
