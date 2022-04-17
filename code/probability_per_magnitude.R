library(ggplot2)
library(viridis)
library(tidyr)
library(scales)
library(here)

rm(list=ls()) # removes all objects in the given environment
wd <- here()
data_dir <- paste0(wd, "/data/")
figure_dir <- paste0(wd, "/figures/")


# index <- 6:23 - 5
# 
# # create a gradent of lambda values
# lam <- seq(0.002, 0.032, 0.00001) # based on Kuo & Ochman (2009), with a slightly smaller lower end
# t <- 4000 # present day diversity, 4000 My of evolution
# df <- data.frame(lam) # initialize dataframe
# auc <- c()


# for (i in index){
#   
#   mag <- i + 5
#   
#   # 10^mag is E(S_t)
#   ep <- function(lam) 1 - (log(10^mag))/(t*lam)
#   
#   # find x intercept 
#   intercept <- optimize(function(lam) 0 - abs(ep(lam)), c(0,0.04), maximum = T)$maximum
#   
#   # integrate curve for specified parameter range
#   area <- integrate(ep, intercept, 0.032)
#   auc <- c(auc, area$value) 
#   # column bind to data frame
#   df <- cbind(df,ep(lam))
#   
#   names(df)[i] <- paste0("ep",mag)
# }
# 
# area_p <- rep(0, length(auc)-1)
# # parameter space area
# for (i in 1:length(area_p)){
#   area_p[i] <- auc[i]-auc[i+1]
# }
# area_prob <- area_p/(sum(area_p))
# 
# area_prob
# mags <- 6:22
# plot(x=mags, y=area_prob)



probabilities <- function(q, p, MEE=FALSE){
  t <- 4000
  lam <- seq(0.002, 0.032, 0.00001)
  df <- data.frame(lam) # initialize dataframe
  
  pal <- viridis_pal()(7)
  auc <- c()
  
  index <- 6:23 - 5
  
  
  for (i in index){
    
    mag <- i + 5
    
    # 10^mag is E(S_t)
    
    if (!MEE) {
      rel_ext <- function(lam) 1 - (log(10^mag))/(t*lam)
      
    }else{
    
      rel_ext <- function(lam) 1 - log(10^mag*(1-p)^-1*(1-q*p)^-5)/(lam*3999)
      
    }
    # ext <- rel_ext(lam, mag, q, p)
    ext <- rel_ext(lam)
    
    intercept <- optimize(function(lam) 0 - abs(rel_ext(lam)), c(0,0.04), maximum = T)$maximum
    
    # integrate curve for specified parameter range
    area <- integrate(rel_ext, intercept, 0.032)
    auc <- c(auc, area$value)
    
    # column bind to data frame
    df <- cbind(df,ext)
    names(df)[i] <- paste0("ep",mag)
  }
  

  
  # parameter space probabilities
  area_p <- rep(0, length(auc)-1)
  # parameter space area
  for (i in 1:length(area_p)){
    area_p[i] <- auc[i]-auc[i+1]
  }
  area_prob <- area_p/(sum(area_p))
  mags <- 6:22
  # plot(x=mags, y=area_prob)
  

  return(area_prob)
}
y0 <- probabilities(0, 1, F)

y1 <- probabilities(1, 0.9, T)
y2 <- probabilities(0.1, 0.9, T)


df <- data.frame(mags = 6:22, E_S=y0) 

p1 <- (ggplot(df, aes(x=mags, y=E_S)) + 
  geom_point(size= 2.0) +
  # geom_smooth(size = 1.0) +
  scale_y_continuous('Probability') +
  scale_x_continuous('Species Richness',
                     breaks=c(6,9, 12, 15, 18, 21),labels = math_format(10^.x)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # legend.position = "none",
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, face = "bold")))

fig3_dir <- paste(figure_dir, "figure3.png", sep = "")
ggsave(plot = p1, filename = fig3_dir, width = 7, height = 5)

# library(drc)
# library(nlme)
# library(aomisc)
# 
# model <- drm(E_S ~ mags, fct = DRC.expoDecay(),
#              data = df)
# 
# model <- lm(E_S ~ poly(mags,2), data=df)
# 
# summary(model)
# 
# plot(x=df$mags, y=df$E_S)
# points(x=df$mags, fitted(model), col='red')

# 
# df <- data.frame(mags = mags, E_S=y0, q1=y1, q_1=y2)
# 
# df1 <- gather(df, key='model', value='probability', -mags)

# ggplot(df1, aes(x=mags, y=probability, group=model)) + 
#   geom_point(aes(shape=model, color=model), size = 2.0) +
#   scale_color_manual(values = c('black', '#6565FE', '#FE6565')) +
#   scale_y_continuous('Probability') +
#   scale_x_continuous('Species Richness',
#                      breaks=c(6,9, 12, 15, 18, 21)) +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         # legend.position = "none",
#         axis.text = element_text(size = 12),
#         axis.line = element_line(colour = "black"),
#         axis.title = element_text(size = 14, face = "bold"))

