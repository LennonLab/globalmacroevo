## CladeRatesMassExtinction.R
## Ford Fishman
# Load packages and functions, set environment
library(scales)
library(MASS)
library(grid)
library(ggplot2)
library(tidyr)
library(png)
source("~/GitHub/MicroSpeciation/code/functions.R")

# Mass extinction events
GOE <- 2450
OrdovicianSilurian <- 445
LateDevonian <- 375
PermianTriassic <- 252
TriassicJurassic <- 201
CretaceousPaleogene <- 66

# vector of all extinction events
EE <- c(GOE, OrdovicianSilurian, LateDevonian, PermianTriassic,TriassicJurassic, CretaceousPaleogene)

time <- 4000
clades <- list(c(1))
lambda <- 0.015
epsilon <- 0.5
S_total <- rep(0, time)
S_total[1] <- 1.0 
pInit <- 0.001 # per clade probability of a clade creating another clade
q <- c(1.0) # susceptibility to mass extinction

# timestep function - discrete time speciation
timestep <- function(clades, t, q, EE){
  
  St_1 <- clades # diversity of clades at previous timestep
  numClades0 <- length(St_1) # number of clades before new clade
  p <- pInit * numClades0 # probability of forming a new clade
  
  if (runif(n=1)<=p){ # clade has new rate
    
    cladeIsDead <- TRUE
    
    while(cladeIsDead){ # is the selected clade extinct
      i <- sample(x = 1:numClades0, size=1)
      cladeIsDead <- St_1[i] < 0
    }
    
    q[numClades0+1] <- runif(n=1) 
    St_1[numClades0+1] <- 1.0
    
  }
  
  numClades1 <- length(lambda) # number of clades
  ep <- rnorm(n = length(numClades0), mean = epsilon, sd = 0.5) # give an epsilon near the actual epsilon for all simulations per time step
  ep <- ifelse(ep<0,0,ep)
  mu <- lambda * ep 
  r <- lambda - mu # diversification/my
  
  St <- St_1*exp(r)*(1-EE*q) # diversity of clades after this timestep
  St <- ifelse(St<1,0,St) # if clade richness is below 1, clade is extinct
  
  return( list(q=q,St=St) )
  
}



for(t in 2:4000) {
  mya <- 4001 - t # how many millions years ago is it?
  
  if (mya %in% EE) {
    ext <- 0.9 # percent of taxa that go extinct from mass extinction
  } else {
    ext <- 0
  }
  
  simList <- timestep(clades=clades[[t-1]], t=t, q=q, EE=ext)
  # update parameters and diversity counts
  q <- simList$q
  clades[[t]] <- simList$St
  S_total[t] <- sum(clades[[t]])
}


totalClades <- max(lengths(clades)) # number of clades with independent rates



l1 <- lapply(clades, 'length<-', max(lengths(clades)))
m1 <- matrix(unlist(l1), ncol = totalClades, nrow = time, byrow = TRUE)
r.S_total <- rev(S_total) # put into mya
df <- as.data.frame(m1)

dominantClades <- unique(colnames(df)[apply(df, 1, which.max)])

df$time <- time:1 -1
df1 <- gather(df, key = "clade", value = "richness", -time)
df1$dominant <- df1$clade %in% dominantClades

p1 <- ggplot(df1, aes(x=time, y=log10(richness), group=clade, color = "Non-Dominant Clade", linetype="Non-Dominant Clade")) +
  geom_line(data = subset(df1, dominant==T),aes(x=time, y=log10(richness), color = "Dominant Clade"),size = 2) +
  scale_x_continuous("Millions of Years Ago", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 10)) +
  scale_y_continuous("Taxon Diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels =  math_format(10^.x), limits = c(0,15)) +
  stat_function(fun=function(time) log10(exp((lambda-lambda*epsilon)*(4000-time))), aes(linetype = "Expected Diversity",color = "Expected Diversity")) +
  geom_vline(xintercept = EE, size = 4, alpha = 0.4) +
  scale_linetype_manual("",
                        limits = c("Dominant Clade","Non-Dominant Clade", "Expected Diversity", "Total Diversity"),
                        breaks = c("Expected Diversity","Total Diversity","Dominant Clade","Non-Dominant Clade"),
                        values = c("solid", "solid", "solid","dotted"), drop = F, guide = 'legend') +
  scale_color_manual("",
                     limits = c("Dominant Clade","Non-Dominant Clade", "Expected Diversity", "Total Diversity"),
                     breaks = c("Expected Diversity","Total Diversity","Dominant Clade","Non-Dominant Clade"),
                     values = c("red","cornflowerblue", "gray18", "black"), drop = F, guide = "legend") +
  annotation_logticks(sides = "lb") +
  guides(color = guide_legend(override.aes = list(size = 0.3))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(5,"pt"),)

p2 <- p1 + geom_line(size = 0.3)
# p2
p3 <- p2 + stat_function(fun = function(time) log10(r.S_total[time+1]), size = 1, aes(linetype = "Total Diversity", color = "Total Diversity")) 
p3


fig3_dir <- paste(figure_dir, "/CladeRatesMassExtinction.png",sep = "")
ggsave(plot = p3, filename = fig3_dir, width = 7, height = 5)
