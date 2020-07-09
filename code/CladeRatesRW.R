## Simulate clade-specific rates by ha
## Ford Fishman

## SETUP ENVIRONMENT
# Load packages
library("png")
library("grid")
library("tidyr")
library("ggplot2")
library("viridis")
library("here")
library("scales")
library("MASS")

rm(list=ls()) # removes all objects in the given environment
wd <- here()
data_dir <- here("data")
figure_dir <- here("figures")



# Set up parameters
time <- 4001
# set up matrix with columns being clade (only 1 to start) and rows representing species over timr
# species <- matrix(S1 = rep(0,time),nrow = time, ncol = length(epsilon)) 

# starting species has the following relative extinction rate:
# species[1, "S1"] = 1
species <- list(c(1))
S_total <- double(time)
S_total[1] <- 1.0 # only 1 starting species
lambda <- c(0.015) # speciation/my, initial value for S1
lams <- list(c(0.015))
mu <- c(0.007)
mus<- list(c(0.007))
pInit <- 0.001 # per clade probability of a clade creating another clade
walk <- 0.0001
# function similating diversification process

  ## Arguments:
  # lambda - vector of speciation rates for all clades for the previous timestep
  # mu - vector of extinction rates for all clades for the previous timestep
  # clades - the richness levels of all clades organized by time and clade
  # ind - the current time

timestep <- function(lambda, mu, species){
  
  St_1 <- species # diversity of clades at previous timestep
  numClades0 <- length(St_1) # number of clades before new clade 
  p <- pInit * numClades0 # probability of forming a new clade
    
  if (runif(n=1)<=p){ # clade has new rate
    
    cladeIsDead <- TRUE
    
    while(cladeIsDead){ # is the selected clade extinct
      i <- sample(x = 1:numClades0, size=1)
      cladeIsDead <- St_1[i] < 0
    }
    
    lambda[numClades0+1] <- lambda[i]
    mu[numClades0+1] <- mu[i]
    St_1[numClades0+1] <- c(1.0) # new clade forms
    
  }
  
  numClades1 <- length(lambda) # number of clades
  testLambda <- runif(numClades1) 
  testMu <- runif(numClades1) 
  
  lambda <- ifelse(testLambda < 0.33, lambda + walk, 
                   ifelse(testLambda < 0.66, lambda, lambda - walk))
  mu <- ifelse(testMu < 0.33, mu + walk, 
                   ifelse(testMu < 0.66, mu, mu - walk))
  lambda <- ifelse(lambda<0,0,lambda)
  lambda <- ifelse(lambda>0.031,0.031,lambda)
  mu <- ifelse(mu<0,0,mu)
  
  r <- lambda - mu # diversification/my

  St <- St_1*r + St_1 # diversity of clades after this timestep
  St <- ifelse(St<1,0,St) # if clade richness is below 1, clade is extinct
   
  
  return(
    list(
      lambda=lambda,
      mu=mu,
      St=St,
      r=r
    )
  )
  
}



# run the function 
for (ind in 2:time){
  
  simList <- timestep(lambda=lambda, mu=mu, species=species[[ind-1]])
  # update parameters and diversity counts
  lambda <- simList$lambda
  mu <- simList$mu
  species[[ind]] <- simList$St
  lams[[ind]] <- lambda
  mus[[ind]] <- mu
  S_total[ind] <- sum(species[[ind]])
}
summary(S_total)

(totalClades <- max(lengths(species)) )# number of clades with independent rates



start <- Sys.time()
l1 <- lapply(species, 'length<-', max(lengths(species)))
m1 <- matrix(unlist(l1), ncol = totalClades, nrow = time, byrow = TRUE)
end <- Sys.time()
end-start

r.S_total <- rev(S_total) # put into mya
df <- as.data.frame(m1)

dominantClades <- unique(colnames(df)[apply(df, 1, which.max)])

# for (t in 1:nrow(df)){
#   a <- df[t,]
#   b <- ifelse(is.na(a),0,df)
#   print(max(a, na.rm = T))
#   
# }

df$time <- time:1 -1
df1 <- gather(df, key = "clade", value = "richness", -time)
df1$dominant <- df1$clade %in% dominantClades

# p1 <- ggplot(df1, aes(x=time, y=log10(richness), group=clade)) +
#   geom_line(size = 0.1, aes(color = "Non-Dominant Clade", linetype="Non-Dominant Clade")) +
#   scale_x_reverse("Time (Mya)", expand = c(0,0), breaks = c(0,1000,2000,3000,4000)) +
#   scale_y_continuous("Taxon Diversity", breaks = c(3,6,9,12), expand = c(0,0), labels =  math_format(10^.x)) +
#   stat_function(fun = function(time) log10(r.S_total[time+1]), size = 1, aes(linetype = "Total Diversity", color = "Total Diversity")) +
#   stat_function(fun=function(time) log10(exp((lams[[1]][1]-mus[[1]][1])*(4000-time))), aes(linetype = "Expected Diversity",color = "Expected Diversity")) +
#   scale_linetype_manual("", breaks = c("Non-Dominant Clade", "Expected Diversity", "Total Diversity"),values = c("solid", "dashed", "dotted")) +
#   scale_color_manual("",breaks = c("Non-Dominant Clade", "Expected Diversity", "Total Diversity"),values = c("cornflowerblue", "gray18", "black")) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14, face = "bold"),
#         axis.ticks = element_line(size = 1),
#         axis.ticks.length = unit(5,"pt"),)
# p1

# only expected value
p1 <- ggplot(df1, aes(x=time, y=log10(richness), group=clade, color = "Non-Dominant Clade", linetype="Non-Dominant Clade")) +
  # geom_line(size = 0.1) +
  geom_line(data = subset(df1, dominant==T),aes(x=time, y=log10(richness), color = "Dominant Clade"),size = 1) +
  scale_x_reverse("Millions of Years Ago", expand = c(0,0), breaks = c(0,1000,2000,3000,4000)) +
  scale_y_continuous("Taxon Diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels =  math_format(10^.x), limits = c(0,15)) +
  # stat_function(fun = function(time) log10(r.S_total[time]), size = 1, aes(linetype = "Total Diversity", color = "Total Diversity")) +
  stat_function(fun=function(time) log10(exp((lams[[1]][1]-mus[[1]][1])*(4000-time))), aes(linetype = "Expected Diversity",color = "Expected Diversity")) +
  # scale_linetype_manual("", labels = c("Expected Diversity","Non-Dominant Clade",  "Total Diversity"),values = c("dashed", "solid",  "dotted"), drop = F) +
  # scale_color_manual("",labels = c("Expected Diversity","Non-Dominant Clade",  "Total Diversity"),values = c("gray18", "cornflowerblue", "black"), drop = F) +
  scale_linetype_manual("",
                        limits = c("Dominant Clade","Non-Dominant Clade", "Expected Diversity", "Total Diversity"),
                        breaks = c("Expected Diversity","Total Diversity","Dominant Clade","Non-Dominant Clade"),
                        values = c("solid", "solid", "solid","dotted"), drop = F, guide = 'legend') +
  scale_color_manual("",
                     limits = c("Dominant Clade","Non-Dominant Clade", "Expected Diversity", "Total Diversity"),
                     breaks = c("Expected Diversity","Total Diversity","Dominant Clade","Non-Dominant Clade"),
                     values = c("red","cornflowerblue", "gray18", "black"), drop = F, guide = "legend") +
  guides(color = guide_legend(override.aes = list(size = 0.3))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)
# p1
# p2 <- p1
# p2$layers[[1]] <- NULL
# p2
p2 <- p1 + geom_line(size = 0.3)
# p2
p3 <- p2 + stat_function(fun = function(time) log10(r.S_total[time+1]), size = 1, aes(linetype = "Total Diversity", color = "Total Diversity")) 
p3

fig1_dir <- paste(figure_dir, "/CladeRatesRW_onlyexpected.png",sep = "")
ggsave(plot = p1, filename = fig1_dir, width = 7, height = 5)

fig2_dir <- paste(figure_dir, "/CladeRatesRW_nototal.png",sep = "")
ggsave(plot = p2, filename = fig2_dir, width = 7, height = 5)

fig3_dir <- paste(figure_dir, "/CladeRatesRW.png",sep = "")
ggsave(plot = p3, filename = fig3_dir, width = 7, height = 5)
