## Simulating Extinction Events Based on Data 
## Ford Fishman

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)
require(ggplot2, quietly = T)
require(viridis, quietly = T)
require(fitdistrplus, quietly = T)
require(goft, quietly = T)
require(actuar, quietly = T)
source("code/functions.R")


## simulations
# extract information about intensitera_lengths and era length from the extinction data
ext_list <- importRohdeMuller()
avg_ext <- ext_list$avg_ext # extinction intensity flattened across duration of era
peak_ext <- ext_list$peak_ext # extinction intensity concentrated on the last year of an era
era_lengths <- ext_list$era_lengths # vector of era lengths
era_means <- ext_list$era_means # vector of extinction levels of eras
original <- ext_list$original # original form of extinction data: not concentrated, not averaged out
EeStart <- as.integer(names(avg_ext[1])) # first year with extinction data

## Simulate simulation events for first 3.5 billion years

# plot data to see distribution
ggplot(NULL,aes(x=era_lengths))+geom_dotplot(size=0.1) 
ggplot(NULL,aes(x=era_means))+geom_dotplot()
lm1 <- lm(era_lengths~era_means) # test for correlation between era length and extinction
summary(lm1) # no significant correlation

# fit extinction to a distribution
fit.gamma <- fitdist(era_means,"gamma") # gamma fit
summary(fit.gamma)
plot(fit.gamma)
fit.gamma.D <- dgamma(1:50,shape = fit.gamma$estimate[1],rate = fit.gamma$estimate[2])

fit.invgauss.means <- fitdist(era_means,"invgauss",start=list(mean = 10, shape = 12)) # inverse gaussian fit
plot(fit.invgauss.means)
hist(era_means,freq = F, breaks = 50)
fit.invgauss.means.D <- dinvgauss(1:50,mean = fit.invgauss.means$estimate[1],shape = fit.invgauss.means$estimate[2])
lines(fit.invgauss.means.D,col = "violet") # a better fit
lines(fit.gamma.D,col = "red")

# fitting era lengths to a distribution
t.lengths <- table(era_lengths)
df.lengths <- as.data.frame(t.lengths)
o_freq.lengths <- c(0,df.lengths$Freq,rep(0,20))
descdist(era_lengths,discrete = T)

fit.nbinom <- fitdist(era_lengths,"nbinom",discrete = T) # negative binomial
fit.nbinom.D <- dnbinom(1:10,size=fit.nbinom$estimate[1],mu=fit.nbinom$estimate[2])

# test with Chi square
# fit.nbinom.D <- dnbinom(0:30,size=fit.nbinom$estimate[1],mu=fit.nbinom$estimate[2])
# chisq.test(o_freq.lengths, p=fit.nbinom.D)

fit.poisson <- fitdist(era_lengths,"pois",discrete = T) # poisson
fit.poisson.D <- dpois(1:10,lambda = fit.poisson$estimate)
plot(fit.poisson)

fit.gamma <- fitdist(era_lengths,"gamma",discrete = T) # gamma 
fit.gamma.D <- dgamma(1:10,shape = fit.gamma$estimate[1],rate = fit.gamma$estimate[2])
plot(fit.gamma)

fit.invgauss <- fitdist(era_lengths,"invgauss", start=list(mean = 2, shape = 4),discrete = T) # inverse gaussian
summary(fit.invgauss)
fit.invgauss.D <- dinvgauss(1:10,mean = fit.invgauss$estimate[1],shape = fit.invgauss$estimate[2])
# Display the various fits
hist(era_lengths,prob=T,breaks=10) 
lines(fit.nbinom.D, col = "blue") 
lines(fit.invgauss.D,col="violet")
lines(fit.poisson.D, col = "red")
lines(fit.gamma.D,col="green")

# The actual simulation 
years = 0
sim_lengths.pool <- round(rinvgauss(3000, mean = fit.invgauss$estimate[1],shape = fit.invgauss$estimate[2])) # simulate era durations with inverse gaussian fit distribution
sim_lengths <- c()

# Make sure the simulated eras last the correct amount of time (~3500 My)
for (i in sim_lengths.pool){
  if (years >= EeStart-1) break
  sim_lengths <- c(sim_lengths, i)
  years=sum(sim_lengths)
}
a <- length(sim_lengths)
sim_lengths[a] <- (EeStart-1) - sum(sim_lengths[1:a-1])
# for each era, simulate a corresponding extinction intensity
sim_intensity <- rinvgauss(sim_lengths,mean = fit.invgauss.means$estimate[1],shape = fit.invgauss.means$estimate[2])

# make sure that every year in the same era is assigned the same extinction intensity
extinction_max <- c()
for(i in 1:a){
  if (sim_lengths[i]!=0){
    extinction_max <- c(extinction_max, rep(sim_intensity[i],sim_lengths[i]))
  }
}

# plot the simulated extinction levels with the actual extinction data
ggplot(NULL, aes(x=4000:1,y=c(extinction_max, original*100))) + 
  geom_line()+
  geom_vline(xintercept = 4000 - EeStart, linetype = "dashed", color = "red") +
  scale_x_reverse("Time (Mya)") +
  ylab("Extinction Intensity (%)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        axis.line.y.right = element_line(color = "grey"), 
        axis.ticks.y.right = element_line(color = "grey")
  )

data <- c(extinction_max, original*100)
presFig.d <- rep(0,4000) # pre-vectorize

for (i in 1:length(data)){
  if (length(data)>(i+20)) {
    j <- i + 20
  }else {
    j <-length(data)
  }
  presFig.d[i] <- mean(data[i:j])
  
}
# plot the simulated extinction levels with the actual extinction data
ggplot(NULL, aes(x=4000:1,y=data)) + 
  geom_line()+
  geom_vline(xintercept = 4000 - EeStart, linetype = "dashed", color = "red") +
  scale_x_reverse("Time (Mya)") +
  ylab("Extinction Intensity (%)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(colour = "black"),
        axis.line.y.right = element_line(color = "grey"), 
        axis.ticks.y.right = element_line(color = "grey")
  )

# plot for presentation
presPlot <- ggplot(NULL, aes(x=4000:1,y=presFig.d)) + 
  geom_line()+
  # geom_vline(xintercept = 4000 - EeStart, linetype = "dashed", color = "red") +
  scale_x_reverse("Time (Mya)") +
  ylab("Extinction Intensity (%)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title=element_text(size=14,face="bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt")
  )

fig_dir1 <- paste(figure_dir, "ExtinctionEventSimPresentation.png", sep = "")
ggsave(plot = presPlot, filename = fig_dir1, width = 7, height = 5)

# test of time correlation of extinction data
summary(extinction_max)
summary(original*100)
times <- EeStart:4000
summary(lm(original~times))
