## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(ggplot2, quietly = T)
require(viridis, quietly = T)
source("code/functions.R")

### VARIABLE SPECIATION AND EXTINCTION 
# Note: reduce number of simulations to reduce run time
# I used 400000 to make the heatmap more smooth

lambda <- runif(n = 400000, min = 0.008, max = 0.031) # speciation rates
epsilon <- runif(n = 400000, min = 0.1, max = 0.9) # mu/lambda
mu <- c() # extinction rates
r <- c() # diversification rates
S <- integer(length = 40000) + 1 # keeps track of species for each set of parameters

for(t in 2:4000) {
  ep <- rnorm(n = epsilon, mean = epsilon, sd = 0.2) # give me an epsilon near the actual epsilon for all simulations per time step
  ep <- ifelse(ep<0,0.01,ep) # if Eb is ever less than 0, make it 0.1
  ep <- ifelse(ep>1,0.9,ep) # upper bound is 0.9
  lam <- rnorm(n = lambda, mean = lambda, sd = 0.01) # do a similar thing for lambda
  lam <- ifelse(lam<0.008,0.008,lam)
  lam <- ifelse(lam>0.031,0.031,lam)
  mu <- lambda * ep # extinction rate
  r <- lambda - mu # diversification rate
  S <- S * exp(r) # discrete exponential growth, pure birth
}

df <- data.frame(lambda, epsilon, mu, r, S)
df <- cbind(df, rounding_scheme(S))

# heatmap 2
(p1 <- ggplot(data = df, aes(x = lambda, y = epsilon, color = roundS)) + 
    geom_point(size = 1.8) +
    xlab(expression(lambda*" (Species/Myr)")) +
    scale_y_continuous(expression(epsilon), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    scale_color_viridis("Log(Richness)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          axis.line = element_blank())
)
# Save figure
fig1_dir <- paste(figure_dir, "VariableSpeciationExtinctionHeatmap.png", sep = "")
ggsave(plot = p1, filename = fig1_dir, width = 7, height = 5)
