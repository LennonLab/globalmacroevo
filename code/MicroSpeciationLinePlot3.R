## Modeling impacts of mass extiction on microbial diversity
## Does not use Sepkoski intensities
## Ford Fishman

library(ggplot2)
library(tidyr)
library(here)
library(stringr)

# Load packages and functions, set environment

setwd(here())
source("code/functions.R")
figure_dir <- "figures/"


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
EE_labels <- c('GOE', 'O-S', 'D', 'P-Tr','Tr-J', 'K-T')
df_EE <- data.frame(year=EE, event=EE_labels)

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

  if (mya %in% EE) {
    
    p <- c(0, 0.1, 0.5, 0.9) # proportion of taxa that go extinct
    
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

df2$intensity <- str_replace(df2$intensity, 'zero', '0%')
df2$intensity <- str_replace(df2$intensity, 'ten', '10%')
df2$intensity <- str_replace(df2$intensity, 'fifty', '50%')
df2$intensity <- str_replace(df2$intensity, 'ninety', '90%')

df3 <- subset(df2, mya>=10)

p1 <- ggplot(df3, aes(x = mya, y = log10(richness), group = intensity)) +  
    geom_line(aes(linetype = intensity)) +
    scale_x_continuous("Millions of Years Ago", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 8)) +
    scale_y_continuous("Taxon Diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels =  math_format(10^.x), limits = c(0,15)) +
    
    scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed","dotdash"),
                        labels = c("0%", "10%","50%", "90%")) +
    annotation_logticks(sides = "b") +
    geom_text(data = subset(df3, mya==10),aes(label = intensity), x = Inf, hjust = 1)+
    geom_text(data=df_EE, aes(x=year,y = 14.5,label = event),size=2, inherit.aes = F)+
    geom_segment(data=df_EE,aes(x = year, xend = year, y = 0, yend=14),size = 4, alpha = 0.4, inherit.aes = F) +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),
        plot.margin = unit(c(1,3,1,1), "lines")
  )
p1

fig1_dir <- paste0(figure_dir, "ExtinctionLinePlot.png")
ggsave(plot = p1, filename = fig1_dir, width = 7, height = 5)
