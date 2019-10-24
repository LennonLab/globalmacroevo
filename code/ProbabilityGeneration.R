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
lambda <- runif(n = sims, min = 0.008, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.1, max = 0.9) # mu/lambda
mu <- lambda * epsilon # extinction rates
r <- lambda - mu # diversification rate
t <- 1
Time <- 4000 # total time to present 
# alpha and beta are terms in the probability generating function
alpha <- (mu*(exp(r*Time)-1))/(lambda*exp(r*Time)-mu)
beta <- alpha*(lambda/mu)

## probability generating function
pgf <- function(S){
  # if only testing the probability of exactly one species richness, 
  # the process is shorter
  if (length(S) == 1){ 
    Prob <- (1-alpha)*(1-beta)*beta^(S-1) # actual PGF from Raup 1985 (A17)
    Prob <- ifelse(Prob<0, 0, Prob) # Probability axioms state they can't be < 0
    RelProb <- Prob/max(Prob) 
    NormProb <- (Prob - mean(Prob))/sd(Prob)
    
  } else { # if testing probability of reaching several different richnesses
    # Vectorize to reduce runtime
    Prob <- double(length(alpha))
    RelProb <- double(length(alpha))
    NormProb <- double(length(alpha))
    
    for (i in 1:length(alpha)){
      prob_set <- (1-alpha[i])*(1-beta[i])*beta[i]^(S-1)
      Prob[i] <- sum(prob_set)
      NormProb[i] <- mean(prob_set*S)
    }
    Prob <- ifelse(Prob<0, 0, Prob)
    RelProb <- Prob/max(Prob)
  }
  return(data.frame(Prob, RelProb, NormProb))
}

## plot relative, normalized, or unaltered probabilities
probPlot <- function(df, type = "RelProb"){
  type <- ifelse(type!="RelProb"&type!="NormProb"&type!="Prob","RelProb",type)
  prob_scale <- ifelse(type=="RelProb","Relative\nProbability",
                  ifelse(type=="NormProb","Normalized\nProbability","Probability"))
  return(ggplot(data = df, aes(x = lambda, y = epsilon, color = get(type))) + 
    geom_point(size = 0.4) +
    xlab(expression(lambda*" (Species/Myr)")) +
    scale_y_continuous(expression(epsilon), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    scale_color_viridis(prob_scale) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          axis.line = element_blank()))
}
# Probability of having exactly 10^12 species
S12 <- 10^12
dfS12 <- data.frame(lambda, epsilon, pgf(S = S12))
pS12 <- probPlot(dfS12, type = "RelProb")


# Probability of being in range for Locey & Lennon 2016
S12_mult <- seq(0.385, 3.43, by = 0.005)*10^12 
dfS12_mult <- data.frame(lambda, epsilon, pgf(S = S12_mult))
pS12_mult <- probPlot(dfS12_mult, type = "NormProb")

# Probability of being in range for Louca et al. 2019
S6_mult <- seq(2.2, 4.3, by = 0.01)*10^6
dfS6_mult <- data.frame(lambda, epsilon, pgf(S = S6_mult))
pS6_mult <- probPlot(dfS6_mult, type = "NormProb")

