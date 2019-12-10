## Ford Fishman
## December 10, 2019

## Construct richness confidence intervals for stochastic birth-death processes
## with different sets of parameter values (speciation, extinction rates).

## SETUP ENVIRONMENT
rm(list=ls())
wd <- "~/Documents/LennonLab/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(R.utils, quietly = T)
require(ggplot2, quietly = T)
require(viridis, quietly = T)

sims <- 100000 # some arbitrary amount of simulations
lambda <- runif(n = sims, min = 0.007, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.0, max = 0.9) # mu/lambda
mu <- lambda * epsilon # extinction rates
r <- lambda - mu # diversification rate
Time <- 4000 # total time to present 
ES <- exp(r*Time) # expected richness
k <- 0 # value to test tail functions against
# subset to constrain complexity
cutoff <- 10^16
lambda <- lambda[ES<cutoff]
epsilon <- epsilon[ES<cutoff]
mu <- mu[ES<cutoff]
r <- r[ES<cutoff]

# from Bailey (1964; eq. 8.46)
beta <- lambda*(exp(r*Time)-1)/(lambda*exp(r*Time) - mu)

# Finds richness confidence intervals for given parameters
CIs <- function(){
  timesup <- F # is time limit up?
  ptm <- proc.time()[3] # start recording time
  while (!timesup&(RoundingLowAll|RoundingHighAll)) { # if time isn't up and CI's aren't all found yet
    # create a random richness within a feasible range for the richness estimates
    k <- runif(n=1,min=0,max=10) * 10^runif(n=1,min=0,max=22) 
    if (RoundingLowAll&RoundingHighAll){ # if some parameter set doesn't have full CI
      # probability that richness at a given time is greater than k
      P_greater <- beta^(k-1)
      # probability that richness at a given time is less than k
      P_lower <- 1 - beta^(k-1)
      # if prob. is 0.975 that process should produce a richness greater than k
      atLowerBound <- round(P_greater, digits = 3)==0.975 
      # if prob. is 0.975 that process should produce a richness less than k
      atUpperBound <- round(P_lower, digits = 3)==0.975
      RoundingLow <- ifelse(atLowerBound|!RoundingLow==F,F,T) # is lower bound found yet
      RoundingHigh <- ifelse(atUpperBound|!RoundingHigh==F,F,T) # is upper bound found yet
      klow <- ifelse(atLowerBound,k,klow) # if yes to above, set lower bound
      khigh <- ifelse(atUpperBound,k,khigh) # same for upper bound
      
      # if only need lower bounds
    } else if(RoundingHighAll){ 
      P_lower <- 1 - beta^(k-1)
      atUpperBound <- round(P_lower, digits = 3)==0.975
      RoundingHigh <- ifelse(atUpperBound|!RoundingHigh==F,F,T)
      khigh <- ifelse(atUpperBound,k,khigh)
      
      # if only need upper bounds
    } else if(RoundingLowAll){
      P_greater <- beta^(k-1)
      atLowerBound <- round(P_greater, digits = 3)==0.975
      RoundingLow <- ifelse(atLowerBound|!RoundingLow==F,F,T)
      klow <- ifelse(atLowerBound,k,klow)
    }
    if (! NA %in% klow) {RoundingLowAll <- F}
    if (! NA %in% khigh) {RoundingHighAll <- F}
    if (proc.time()[3]-ptm>60*5) {timesup <- T} # don't run process for more than 10s 
  }

  return(list(klow=klow, khigh=khigh)) # return bounds as a list of 2 vectors
}

newSims <- length(r) # after controlling for numbers that are too large, how many simulations are run
klow <- rep(NA, newSims) # pre-vectorize with NAs
khigh <- rep(NA, newSims)

RoundingLow = !logical(newSims) # has the lower boundary been found? a vector of booleans (TRUE)
RoundingHigh = !logical(newSims) # has the upper boundary been found?
RoundingLowAll = TRUE 
RoundingHighAll = TRUE

ptm <- proc.time() # start time
klist <- CIs() # run main process
proc.time() - ptm # check process run time#
klows <- klist[[1]]
khighs <- klist[[2]]
newSims
# summary(klows)
# summary(khighs)
# hasData <- !is.na(klows)&!is.na(khighs)
# klows <- klows[hasData]

# create df for plot
df <- data.frame(lambda, epsilon, mu, r, klows, khighs)

# Locey & Lennon (2016) interval
low12 <- 0.385*10^12
high12 <- 3.43*10^12
# Louca et al. 2019 interval
low6 <- 2.2*10^6
high6 <- 4.3*10^6

# is there overlap between the CIs and above intervals?
df$CI12 <- !(klows<low12 & khighs<low12) & !(klows>high12 & khighs>high12) 
df$CI6 <- !(klows<low6 & khighs<low6) & !(klows>high6 & khighs>high6) 

df1 <- subset(df, !is.na(CI6)) # remove NAs
df1$CI <- ifelse(df1$CI12,"CI12",ifelse(df1$CI6,"CI6","none")) # categorize simulations

# plot figure
p1 <- ggplot(data = df1, aes(x = lambda, y = epsilon, color = CI)) + 
  geom_point(size = 1.4) +
  xlab(expression(lambda*" (Species/Myr)")) +
  scale_y_continuous(expression(epsilon), breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  scale_color_manual(breaks = c("CI6","CI12","none"),
                     values = c("dodgerblue3","mediumseagreen","darkslategrey")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        axis.line = element_blank())
# Save figure
fig_dir <- paste(figure_dir, "ConfidenceIntervals.png", sep = "")
ggsave(plot = p1, filename = fig_dir, width = 7, height = 5)
