## Ford Fishman
## December 10, 2019

## Use a probability mass function to give probability of seeing certain levels of microbial 
## richness giving a set of parameters.
#
## Assumes stochastic birth-death process starting with 1 species.
#
## Traditionally, probabilities from a PMF are summed over a desired range of discrete outcomes to 
## calculate the probability of outcomes from that range occuring. However, given the scale of the
## ranges we are considering, vector lengths exceed 10^12, so the integrate() function is used in a
## similar fashion to how one would approach probability calculation for a proability density function
## (PDF). Testing showed that even at the scales of 10^6, integration proved to give similar outcomes 
## to summation. 

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/Documents/LennonLab/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(ggplot2, quietly = T)
require(viridis, quietly = T)

## define parameters alpha and beta
sims = 400000 # number of simulations
lambda <- runif(n = sims, min = 0.001, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.0, max = 0.9) # mu/lambda
mu <- lambda * epsilon # extinction rates
r <- lambda - mu # diversification rate
Time <- 4000 # total time to present

# alpha and beta are terms in the probability generating function
alpha <- (mu*(exp(r*Time)-1))/(lambda*exp(r*Time)-mu) # probability of extinction 
beta <- alpha*(lambda/mu) # modifier of alpha


## probability mass function
pmf_integration <- function(Smin, Smax){
  # make vector in advanced to reduce runtime
  prob <- double(length(lambda)) 
  # need to use for loop, integrate() doesn't like vectors
  for (i in 1:length(lambda)){ 
    # This pmf is taken from Raup (1985; Eq. A17), which is in turn from Bailey (1964; eq 8.46)
    pmf <- function(S) (1-alpha[i])*(1-beta[i])*beta[i]^(S-1) 
    # integrating overall interval of given richness estimation
    prob[i] <- integrate(pmf, Smin, Smax)$value
  }
  return(prob)
}

## plot probabilities 
probPlot <- function(df){
  return(ggplot(data = df, aes(x = lambda, y = epsilon, color = prob)) + 
    geom_point(size = 0.4) +
    xlab(expression(lambda*" (Species/Myr)")) +
    scale_y_continuous(expression(epsilon), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    scale_color_viridis("Probability") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          axis.line = element_blank()))
}
# Probability of having exactly 10^12 species
dfS12 <- data.frame(lambda, epsilon, prob = pmf_integration(0.385*10^12, 3.43*10^12))
(pS12 <- probPlot(dfS12))
# Save figure
figS12_dir <- paste(figure_dir, "probS12_stochastic.png", sep = "")
ggsave(plot = pS12, filename = figS12_dir, width = 7, height = 5)

# Probability of being in range for Louca et al. 2019
dfS6 <- data.frame(lambda, epsilon, prob = pmf_integration(2.2*10^6, 4.3*10^6))
(pS6 <- probPlot(dfS6))
figS6_dir <- paste(figure_dir, "probS6_stochastic.png", sep = "")
ggsave(plot = pS6, filename = figS6_dir, width = 7, height = 5)

# Probabilities for both sets of estimates on same plot
df6_12 <- data.frame(lambda, epsilon, prob = dfS12$prob + dfS6$prob)
pS6_12 <- probPlot(df6_12)
figS6_12_dir <- paste(figure_dir, "probS6_12_stochastic.png", sep = "")
ggsave(plot = pS6_12, filename = figS6_12_dir, width = 7, height = 5)

## variance from Raup 1985 A22, via Bailey (1964; eq 8.49)
variance <- (lambda + mu)/(lambda - mu) * exp((lambda - mu)*Time)*(exp((lambda - mu)*Time) - 1)
expectS <- exp(r * Time)
dfvar <- data.frame(lambda, epsilon, cov = (variance)^0.5/expectS)
(pvar <- ggplot(data = dfvar, aes(x = lambda, y = epsilon, color = cov)) + 
    geom_point(size = 0.4) +
    xlab(expression(lambda*" (Species/Myr)")) +
    scale_y_continuous(expression(epsilon), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    scale_color_viridis("Coefficient\nof Variation") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          axis.line = element_blank()))

# Save figure
var_fig_dir <- paste(figure_dir, "CoV_stochastic.png", sep = "")
ggsave(plot = pvar, filename = var_fig_dir, width = 7, height = 5)

## TESTING INTEGRATION VS SUMMING
## it seems that for the scales we are dealing with, the probabilities aren't different enough
## to worry that PMF's should be summed over integers and not integrated

# lambda <- 0.05
# epsilon <- 0.1
# mu <- lambda*epsilon
# r <- lambda - mu
# Time <- 400
# log10(exp(r*Time)) # check general scale of number
# alpha <- (mu*(exp(r*Time)-1))/(lambda*exp(r*Time)-mu) # probability of extinction 
# beta <- alpha*(lambda/mu) 
# func <- function(S) (1-alpha)*(1-beta)*beta^(S-1) #pmf
# Srange <- (0.385*10^6):(3.43*10^7) # a appropriate for the scale in question
# integrate(func, 0.385*10^6, 3.43*10^7)$value # the integration-based probability
# sum(func(Srange)) # the sum-based probability

