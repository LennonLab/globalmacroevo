## Expectations_heatmap.R
## Ford Fishman
# Load packages
library(png)
library(grid)
library(kableExtra)
library(tidyr)
library(ggplot2)
library(scales)
library(viridis)
library(here)
# rm(list=ls()) # removes all objects in the given environment
wd <- here()
data_dir <- paste0(wd, "/data/")
figure_dir <- paste0(wd, "/figures/")

# create a gradent of lambda values
lam <- seq(0.004, 0.031, 0.00001) # based on Kuo & Ochman (2009), with a slightly smaller lower end
t <- 4000 # present day diversity, 4000 My of evolution
df <- data.frame(lam) # initialize dataframe

# for S = {10^6, 10^9, 10^12, 10^15, 10^18}
for (i in 2:6){
  
  mag <- i*3
  
  # 10^mag is E(S_t)
  ep <- 1 - (log(10^mag))/(t*lam)
  
  # column bind to data frame
  df <- cbind(df,ep)
  names(df)[i] <- paste0("ep",mag)
}

# have to do 10^23 manually
ep23 <- 1 - (log(10^23))/(t*lam)
df <- cbind(df,ep23)

# for correct positions of labels
y1 <- function(x) -30*x +1 # equation of diagonal label line
xs <- rep(0, 5)
ys <- xs
labels <- xs

ind <- 0
for (i in seq(6,18,3)){
  ind <- ind + 1
  
  # want a different label offset for 10^18
  if (i < 18){
    y2 <- function(x) 1 - (log(10^(i+1.5)))/(4000*x)
    
  } else {
    y2 <- function(x) 1 - (log(10^(i+2)))/(4000*x)
  }
  xx <- optimize(function(x) abs(y1(x)-y2(x)), c(.01,.03))$minimum
  ys[ind] <- y1(xx)
  xs[ind] <- xx
  labels[ind] <- paste0('10^',i)
}

# color palette
pal <- viridis_pal()(7)

# plot
(p1 <- ggplot(df) +
  geom_area(aes(x=lam,y=ep6), fill=pal[2])+
  geom_area(aes(x=lam,y=ep9), fill=pal[3])+
  geom_area(aes(x=lam,y=ep12), fill=pal[4])+
  geom_area(aes(x=lam,y=ep15), fill=pal[5])+
  geom_area(aes(x=lam,y=ep18), fill=pal[6])+
  geom_area(aes(x=lam,y=ep23), fill=pal[7])+
  annotate(geom="text", label='10^5',x=0.008, y = y1(0.008), color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label=labels,x=xs, y = ys, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^23',x=0.023, y = y1(0.023), color = 'black', size = 6, parse=TRUE) +
  scale_x_continuous(expression("Speciation Rate, "*lambda*" (Species/Myr)"),expand = c(0,0)) +
  scale_y_continuous("Relative Extinction Rate, \u03B5", 
                     limits = c(0.0, 0.9),
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                     expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill=pal[1]),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.text.align = 1,
        axis.title = element_text(size = 13)))
# Save figures
fig1_dir <- paste(figure_dir, "ExpectationBD_Heatmap.png", sep = "")
ggsave(plot = p1, filename = fig1_dir, width = 7, height = 5)
