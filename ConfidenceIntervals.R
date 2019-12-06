# Confidence intervals
rm(list=ls())
wd <- "~/Documents/LennonLab/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(R.utils, quietly = T)
require(ggplot2, quietly = T)
require(viridis, quietly = T)

sims <- 100000
lambda <- runif(n = sims, min = 0.007, max = 0.031) # speciation rates
epsilon <- runif(n = sims, min = 0.0, max = 0.9) # mu/lambda
mu <- lambda * epsilon # extinction rates
r <- lambda - mu # diversification rate
Time <- 4000 # total time to present 
ES <- exp(r*Time) # expected richness

# subset to constrain complexity
cutoff <- 10^16
lambda <- lambda[ES<cutoff]
epsilon <- epsilon[ES<cutoff]
mu <- mu[ES<cutoff]
r <- r[ES<cutoff]

k <- 0
# from Magallon et al 2002
beta <- (exp(r*Time)-1)/(exp(r*Time) - mu/lambda)

CIs <- function(){
  timesup <- F
  ptm <- proc.time()[3]
  while (!timesup&(RoundingLowAll|RoundingHighAll)) {
    k <- runif(n=1,min=0,max=10) * 10^runif(n=1,min=0,max=22)

    if (RoundingLowAll&RoundingHighAll){
      # probability that richness at a given time is greater than k
      P_greater <- beta^(k-1)
      # probability that richness at a given time is less than k
      P_lower <- 1 - beta^(k-1)
      atLowerBound <- round(P_greater, digits = 3)==0.975
      atUpperBound <- round(P_lower, digits = 3)==0.975
      RoundingLow <- ifelse(atLowerBound|!RoundingLow==F,F,T)
      RoundingHigh <- ifelse(atUpperBound|!RoundingHigh==F,F,T)
      klow <- ifelse(atLowerBound,k,klow)
      khigh <- ifelse(atUpperBound,k,khigh)
    } else if(RoundingHighAll){
      P_lower <- 1 - beta^(k-1)
      atUpperBound <- round(P_lower, digits = 3)==0.975
      RoundingHigh <- ifelse(atUpperBound|!RoundingHigh==F,F,T)
      khigh <- ifelse(atUpperBound,k,khigh)
    } else if(RoundingLowAll){
      P_greater <- beta^(k-1)
      atLowerBound <- round(P_greater, digits = 3)==0.975
      RoundingLow <- ifelse(atLowerBound|!RoundingLow==F,F,T)
      klow <- ifelse(atLowerBound,k,klow)
    }
    if (! NA %in% klow) {RoundingLowAll <- F}
    if (! NA %in% khigh) {RoundingHighAll <- F}
    if (proc.time()[3]-ptm>60*5) {timesup <- T} # dont run process for more than 10s 
  }
  # kvals <- list(klow=klow, khigh=khigh)
  return(list(klow=klow, khigh=khigh))
}
# data_m <- matrix(0,nrow=sims,ncol=5)
newSims <- length(r)
klow <- rep(NA, newSims)
khigh <- rep(NA, newSims)

RoundingLow = !logical(newSims) # has the lower boundary been found? a vector of booleans (TRUE)
RoundingHigh = !logical(newSims) # has the upper boundary been found?
RoundingLowAll = TRUE
RoundingHighAll = TRUE

ptm <- proc.time()
klist <- CIs()
proc.time() - ptm
klows <- klist[[1]]
khighs <- klist[[2]]
newSims
# summary(klows)
# summary(khighs)
# hasData <- !is.na(klows)&!is.na(khighs)
# klows <- klows[hasData]


df <- data.frame(lambda, epsilon, mu, r, klows, khighs)

low12 <- 0.385*10^12
high12 <- 3.43*10^12
# low12 <- 0.00001*(10^12)
# high12 <- 1000*(10^12)
low6 <- 2.2*10^6
high6 <- 4.3*10^6

df$CI12 <- !(klows<low12 & khighs<low12) & !(klows>high12 & khighs>high12) 
df$CI6 <- !(klows<low6 & khighs<low6) & !(klows>high6 & khighs>high6) 

df1 <- subset(df, !is.na(CI6))
df1$CI <- ifelse(df1$CI12,"CI12",ifelse(df1$CI6,"CI6","none"))
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
