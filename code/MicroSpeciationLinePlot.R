## Plotting Microbial Speciation against Extinction Rates
## Ford Fishman
library(scales)
library(MASS)
library(grid)
library(png)
# Load packages and functions, set environment
source("~/GitHub/MicroSpeciation/code/functions.R")

## simulations
ext_list <- importRohdeMuller()
avg_ext <- ext_list$avg_ext
peak_ext <- ext_list$peak_ext
original <- ext_list$original
# timestep function - discrete time speciation
timestep <- function(S, r, ep = 0){
  return(S * exp(r)* (1 - ep))
}

lambda <- 0.020 # speciation
epsilon <- c(0.5, 0.5, 0.5, 0.5) # some extinction rates
S1 <- c(1) # initial species levels 4Ga - first parameter set
S2 <- c(1) # initial species levels 4Ga - second parameter set
S3 <- c(1) # initial species levels 4Ga - third parameter set
S4 <- c(1) # initial species levels 4Ga - fourth parameter set
EeStart <- as.integer(names(avg_ext[1])) # first year with extinction data
EE <- c(rep(0,EeStart-1),original)
# Simulate, but replace random background extinction levels 
# with "known" extinction rates for macroorganisms for "recent" years
for(t in 2:4000) {
  ep <- rnorm(n = epsilon, mean = epsilon, sd = 0.1) # give me an epsilon near the actual epsilon for all simulations per time step
  ep <- ifelse(ep<0,0.01,ep) # if epsilon is ever less than 0, make it 0.1
  ep <- ifelse(ep>1,0.9,ep) # don't want epsilon above 0.9
  mu <- lambda * ep 
  r <- lambda - mu # diversification rate
  
  S1[t] <- timestep(S = S1[t-1], r = r[1]) # discrete exponential growth, pure birth
  S2[t] <- timestep(S = S2[t-1], r = r[2])
  S3[t] <- timestep(S = S3[t-1], r = r[3])
  S4[t] <- timestep(S = S4[t-1], r = r[4])
  if (t >= EeStart){ # when extinction data starts, use it in addition to background extinction
    # ep <- avg_ext[t - EeStart + 1]
    if (t < 4000){
      # ep <- avg_ext[t - EeStart + 1]
      ep <- peak_ext[t - EeStart + 1]
    } else {
      # ep <- median(avg_ext)
      ep <- median(peak_ext)
    }
    ep <- ep*c(0,0.1,0.5,0.9)
    S1[t] <- timestep(S = S1[t-1], r = r[1], ep = ep[1]) # ep % of population dies
    S2[t] <- timestep(S = S2[t-1], r = r[2], ep = ep[2])
    S3[t] <- timestep(S = S3[t-1], r = r[3], ep = ep[3])
    S4[t] <- timestep(S = S4[t-1], r = r[4], ep = ep[4])
  }
}
df <- data.frame(logrichness = c(log10(S1), log10(S2), log10(S3),log10(S4)), time = rep(4000:1, 4), 
                 lambda = c(rep(lambda[1],4000*4)),
                 epsilon = c(rep(epsilon[1],4000), rep(epsilon[2],4000), rep(epsilon[3],4000),rep(epsilon[4],4000)),
                 extinction = rep(EE, 4),
                 relative = c(rep("0%",4000),rep("10%",4000),rep("50%",4000),rep("90%",4000))
                 )
(p1 <- ggplot(data = df, aes(y = logrichness, x = time, group = as.character(relative))) + 
  geom_line(aes(linetype = as.character(relative),color = "Richness")) +
  geom_line(aes(y=extinction*30,color = "Extinction")) +
  scale_x_reverse("Time (Mya)",limits = c(800,0)) +
  scale_y_continuous("Log(Richness)", breaks = seq(0,32,4), expand = c(0,0),
                     sec.axis = sec_axis(~./30, name = "Extinction Event Intensity",labels = percent) ) +
  scale_linetype_discrete() +
  scale_color_manual("",values = c("grey","black")) +
  geom_vline(xintercept = 4000 - EeStart, linetype = "dotted", color = "red") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        axis.line.y.right = element_line(color = "grey"), 
        axis.ticks.y.right = element_line(color = "grey")
))

fig_dir1 <- paste(figure_dir, "ExtinctionEvents_lambdalineplots.png", sep = "")
ggsave(plot = p1, filename = fig_dir1, width = 7, height = 5)

p2 <- ggplot(data = df, aes(y = logrichness, x = time, group = as.character(relative)), color = "black") + 
    geom_line(aes(linetype = as.character(relative))) +
    # geom_line(aes(y=extinction*20,color = "Extinction")) +
    scale_x_reverse("Millions of years ago",limits = c(800,0), expand = c(0,0)) +
    scale_y_continuous("Taxon richness", breaks = seq(0,32,4),minor_breaks = c(seq(0,32,2)), limits = c(8,18),labels =  math_format(10^.x),expand = c(0,0)) +
    # scale_y_continuous("Taxon richness", breaks = seq(4,32,4), labels =  math_format(10^.x),expand = c(0,0),
    #                    sec.axis = sec_axis(~./20, name = "Extinction Intensity\n(Macro-organisms)",labels = percent, breaks = c(0.1,0.3,0.5)) ) +
    scale_linetype_discrete(expression(lambda)) +
    # scale_color_manual("",values = c("dodgerblue4","black")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(5,"pt")
    )
p3 <- ggplot(data = df, aes(y = extinction, x = time, group = as.character(relative))) +
  geom_line(color = "dodgerblue4") +
  scale_x_reverse("Millions of years ago",limits = c(800,0), expand = c(0,0)) +
  scale_y_continuous("Extinction Intensity\n(Macro-organisms)",labels = percent, breaks = c(0.1,0.3,0.5)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt")
  )
fig_dir2 <- paste(figure_dir, "Pres_ExtinctionEvents_VariableImpact.png", sep = "")
# open.
# png(fig_dir2,width = 1400, height = 1000)
# grid.newpage()
p <- rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last")
# grid.draw(rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last"))
png(fig_dir2, width = 7, height = 5, units = "in", res = 1000, pointsize = 12)
# grid.draw(p)
# # grid.(rbind(ggplotGrob(p2), ggplotGrob(p3), size = "last"))
# 
dev.off()

ggsave(filename = fig_dir2, plot = p, width = 7, height = 5)
