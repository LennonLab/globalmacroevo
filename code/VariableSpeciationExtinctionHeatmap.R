## SETUP ENVIRONMENT
require(ggplot2, quietly = T)
require(viridis, quietly = T)
library("here")
# rm(list=ls()) # removes all objects in the given environment
# wd <- "~/GitHub/MicroSpeciation"
data_dir <- here("data")
figure_dir <- here("figures")
# data_dir <- paste(wd, "/data/", sep = "")
# figure_dir <- paste(wd, "/figures/", sep = "")
# getwd()
# setwd(wd)

source(here("code","functions.R"))

### VARIABLE SPECIATION AND EXTINCTION 
# Note: reduce number of simulations to reduce run time
# I used 400000 to make the heatmap more smooth
sims <- 400000
lambda <- runif(n = sims, min = 0.008, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.1, max = 0.9) # mu/lambda
ep <- epsilon 
mu <- c() # extinction rates
r <- c() # diversification rates
S <- integer(length = sims) + 1 # keeps track of species for each set of parameters
#
# pop <- c(1)

for(t in 2:4000) {
  prob <- runif(n=sims, 0, 1) # random variable
  ep <- ifelse(prob<0.33, ep - 0.005, 
               ifelse(prob<0.66, ep, ep + 0.005))
  ep <- ifelse(ep<0,0,ep) # if Eb is ever less than 0, make it 0.1
  ep <- ifelse(ep>0.9,0.9,ep) # upper bound is 0.9
  lam <- rnorm(n = sims, mean = lambda, sd = 0.01) # do a similar thing for lambda
  lam <- ifelse(lam<0.008,0.008,lam)
  lam <- ifelse(lam>0.031,0.031,lam)
  mu <- lambda * ep # extinction rate
  r <- lambda - mu # diversification rate
  S <- S * exp(r) # discrete exponential growth, pure birth
  # pop[t] <- pop[t-1] * exp(r)
}

df <- data.frame(lambda, epsilon, mu, r, S)
df <- cbind(df, rounding_scheme(S))

# heatmap 2
(p1 <- ggplot(data = df, aes(x = lambda, y = epsilon, color = roundS)) + 
    geom_point(size = 1.8) +
    scale_x_continuous(expression(lambda*" (Species/Myr)"),expand = c(0,0)) +
    scale_y_continuous(expression(epsilon), 
                       limits = c(0.1, 0.9),
                       breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                       expand = c(0,0)) +
    scale_color_viridis("Log(Richness)") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(size = 11),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(5,"pt"),
          legend.position = "right",
          axis.line = element_blank(),
    )
)
# Save figure
fig1_dir <- paste(figure_dir, "VariableSpeciationExtinctionAutoHeatmap.png", sep = "")
ggsave(plot = p1, filename = fig1_dir, width = 7, height = 5)


ggplot(NULL, aes(x = 1:4000, y = log10(pop))) + geom_line()
