library(ggplot2)
library(viridis)
library(tidyr)
library(scales)
library(here)

rm(list=ls()) # removes all objects in the given environment
wd <- here()
data_dir <- paste0(wd, "/data/")
figure_dir <- paste0(wd, "/figures/")

probabilities <- function(q, p, MEE=FALSE){
  t <- 4000
  lam <- seq(0.004, 0.030, 0.00001)
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
    area <- integrate(rel_ext, intercept, 0.030)
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
  

  return(area_prob)
}
y0 <- probabilities(0, 1, F)

y1 <- probabilities(1, 0.9, T)
y2 <- probabilities(0.1, 0.9, T)


df <- data.frame(mags = 6:22, E_S=y0, q_1=y1, q_.1=y2) 
df <- gather(df, key='model', value='probability', -mags)

legend_labels <- c("BDE","ME q = 1","ME q = 0.1")

(p1 <- ggplot(df, aes(x=mags, y=probability, group=model)) + 
  geom_point(aes(shape=model, color=model), size= 2.0) +
  scale_y_continuous('Probability') +
  scale_x_continuous('Species diversity',
                     breaks=c(6,9, 12, 15, 18, 21),labels = math_format(10^.x)) +
  scale_color_manual("Model",values = c('#1a7ce5', '#82e619', '#eb1488'), labels = legend_labels) +
  scale_shape_manual("Model",values = c(19, 15, 17), labels = legend_labels) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, face = "bold")))

figS1_dir <- paste(figure_dir, "figureS1.png", sep = "")
ggsave(plot = p1, filename = figS1_dir, width = 7, height = 5)

##### compare ratio of feasible area to infeasible

ratios <- function(MEE=F, p=0.9, q=1.0){
  t <- 4000
  lam <- seq(0.004, 0.030, 0.00001)
  df <- data.frame(lam) # initialize dataframe
  
  ceiling <- function(lam) 1
  
  area <- integrate(Vectorize(ceiling), 0.004, 0.030)
  
  auc <- c(area$value) # area under the first curve is 1, bc the first curve is epsilon = 1 (edge of feasible space)
  
  mags <- c(6, 23)
  
  for (mag in mags){
  
    if (!MEE) {
      rel_ext <- function(lam) 1 - (log(10^mag))/(t*lam)
      
    }else{
      
      rel_ext <- function(lam) 1 - log(10^mag*(1-p)^-1*(1-q*p)^-5)/(lam*3999)
      
    }
    
    ext <- rel_ext(lam)
    
    intercept <- optimize(function(lam) 0 - abs(rel_ext(lam)), c(0,0.04), maximum = T)$maximum
    
    # integrate curve for specified parameter range
    area <- integrate(rel_ext, intercept, 0.030)
    auc <- c(auc, area$value)
    
  }
  
  auc <- c(auc, 0) # add in a zero for the math
  
  area_p <- rep(0, length(auc)-1)
  
  for (i in 1:length(area_p)){
    area_p[i] <- auc[i]-auc[i+1]
  }
  area_prob <- area_p/(sum(area_p))
  names(area_prob) <- c('below 10^6', 'feasible', 'above 10^23')
  return(area_prob)
}

print('BDE:')
ratios(MEE = F)
print('MEE (q=1):')
ratios(MEE = T, q=1)
print('MEE:(q=0.1)')
ratios(MEE = T, q=0.1)
