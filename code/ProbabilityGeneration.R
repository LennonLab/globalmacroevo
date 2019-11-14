## Ford Fishman
## October 22, 2019

## Use a probability generating function to give probability of seeing a certain microbial
## richness giving a set of parameters.
## Assumes stochastic birth-death process starting with 1 species. 

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/Documents/LennonLab/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(ggplot2)
require(viridis)

## define parameters alpha and beta
sims = 400000 # number of simulations
lambda <- runif(n = sims, min = 0.001, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.1, max = 0.9) # mu/lambda
mu <- lambda * epsilon # extinction rates
r <- lambda - mu # diversification rate
t <- 1
Time <- 4000 # total time to present 
# alpha and beta are terms in the probability generating function
alpha <- (mu*(exp(r*Time)-1))/(lambda*exp(r*Time)-mu)
beta <- alpha*(lambda/mu)

## density function
# func <- function(S) (1-alpha)*(1-beta)*beta^(S-1)
## probability generating function
pgf <- function(Smin, Smax){
  prob <- double(length(lambda))
  for (i in 1:length(lambda)){
    func <- function(S) (1-alpha[i])*(1-beta[i])*beta[i]^(S-1)
    prob[i] <- integrate(func, Smin, Smax)$value
  }
  return(prob)
}

## plot relative, normalized, or unaltered probabilities
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
dfS12 <- data.frame(lambda, epsilon, prob = pgf(0.385*10^12, 3.43*10^12))
(pS12 <- probPlot(dfS12))


# Probability of being in range for Louca et al. 2019
dfS6 <- data.frame(lambda, epsilon, prob = pgf(2.2*10^6, 4.3*10^6))
(pS6 <- probPlot(dfS6))

df6_12 <- data.frame(lambda, epsilon, prob = dfS12$prob + dfS6$prob)
probPlot(df6_12)

## variance
variance <- (lambda + mu)/(lambda - mu) * exp((lambda - mu)*Time)*(exp((lambda - mu)*Time) - 1)
expectS <- exp(r * Time)
dfvar <- data.frame(lambda, epsilon, prob = (variance)^0.5/expectS)
(pvar <- probPlot(dfvar))

# Confidence intervals
# lambda <- 0.015
# mu <- 0.007
r <- lambda - mu
beta <- (exp(r*Time)-1)/(exp(r*Time) - mu/lambda)
P_greater <- function(k) {
  return(beta^(k-1))
}
P_lower <- function(k) {
  return(1 - beta^(k-1))
}

stillRoundingLow = TRUE
stillRoundingHigh = TRUE
klow <- 0
khigh <- 10^100
ks <- double(sims)
CIs <- function()
while (stillRoundingLow|stillRoundingHigh) {
  k <- runif(n = 1,min = 1, max = 10) * 10 ^runif(n=1, 5, 23)
  if (stillRoundingLow&stillRoundingHigh){
    lower <- P_greater(k)
    upper <- P_lower(k)
    stillRoundingLow <- ifelse(round(lower, digits = 3)==0.975,FALSE,TRUE)
    stillRoundingHigh <- ifelse(round(upper, digits = 3)==0.975,FALSE,TRUE)
    klow <- ifelse(!stillRoundingLow,k,klow)
    khigh <- ifelse(!stillRoundingHigh,k,khigh)
  } else if(stillRoundingHigh){
    upper <- P_lower(k)
    stillRoundingHigh <- ifelse(round(upper, digits = 3)==0.975,FALSE,TRUE)
    khigh <- ifelse(!stillRoundingHigh,k,khigh)
  } else{
    lower <- P_greater(k)
    stillRoundingLow <- ifelse(round(lower, digits = 3)==0.975,FALSE,TRUE)
    klow <- ifelse(!stillRoundingLow,k,klow)
  }
}
klow
khigh
