## bd_expectations.R
## Ford Fishman

## Generate panels for figure 1 and combine

# Load packages
library(png)
library(grid)
library(kableExtra)
library(tidyr)
library(ggplot2)
library(scales)
library(viridis)
library(here)
library(cowplot)

# Directory information
wd <- here() # project directory
figure_dir <- here("figures")

######## panel a #############
# Levels of diversity across speciation and extinction rates

# create a gradient of lambda values
lam <- seq(0.004, 0.030, 0.00001) # based on Kuo & Ochman (2009), with a slightly smaller lower end
t <- 4000 # present day diversity, 4000 My of evolution
dfA <- data.frame(lam) # initialize dataframe

# create a contour of Eq. 3 for various levels of present-day S
# for S = {10^6, 10^9, 10^12, 10^15, 10^18}
for (i in 2:7){
  
  mag <- i*3 
  
  # 10^mag is E(S_t)
  ep <- function(lam) 1 - (log(10^mag))/(t*lam)
  # add to data frame
  dfA <- cbind(dfA,ep(lam)) 
  names(dfA)[i] <- paste0("ep",mag) 
}

# have to do 10^23 manually
ep23 <- function(lam) 1 - (log(10^23))/(t*lam)
dfA <- cbind(dfA,ep23(lam))
names(dfA)[ncol(dfA)] <- paste0("ep","23")

# for correct positions of labels, place them on a line in their respective bins
y1 <- function(x) -30*x +1 # equation of diagonal label line
xs <- rep(0, 5)
ys <- xs
labels <- xs

ind <- 0
for (i in seq(6,21,3)){
  ind <- ind + 1
  
  # want a different label offset for 10^21
  if (i < 21){
    y2 <- function(x) 1 - (log(10^(i+1.5)))/(4000*x)
    
  } else {
    y2 <- function(x) 1 - (log(10^(i+1)))/(4000*x)
  }
  xx <- optimize(function(x) abs(y1(x)-y2(x)), c(.01,.03))$minimum
  ys[ind] <- y1(xx)
  xs[ind] <- xx
  labels[ind] <- paste0('10^',i)
}

# color palette
pal <- viridis_pal()(8)

# plot
pA <- ggplot(dfA) +
  geom_area(aes(x=lam,y=ep6), fill=pal[2])+
  geom_area(aes(x=lam,y=ep9), fill=pal[3])+
  geom_area(aes(x=lam,y=ep12), fill=pal[4])+
  geom_area(aes(x=lam,y=ep15), fill=pal[5])+
  geom_area(aes(x=lam,y=ep18), fill=pal[6])+
  geom_area(aes(x=lam,y=ep21), fill=pal[7])+
  geom_area(aes(x=lam,y=ep23), fill=pal[8])+
  annotate(geom="text", label='10^5',x=0.008, y = y1(0.008), color = 'white', size = 5.5, parse=TRUE) +
  annotate(geom="text", label=labels,x=xs, y = ys, color = 'white', size = 5.5, parse=TRUE) +
  annotate(geom="text", label='10^23',x=0.023, y = y1(0.023), color = 'black', size = 5.5, parse=TRUE) +
  scale_x_continuous(expression(bold("Speciation rate, "*lambda* " "(My^-1))),
                     limits=c(0.004, 0.030),
                     breaks=seq(0.005, 0.03, 0.005),
                     expand = c(0,0)) +
  scale_y_continuous("Relative extinction rate, \u03B5", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                     expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill=pal[1]),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.text.align = 1,
        axis.title = element_text(size = 16, face = "bold"),
        plot.margin = unit(c(1,2,1,1), "lines"))


######## panel b #############
# probability of diversity outcomes across speciation rates

S <- c(6,9,12,15,18,21,23) # log10 S values

# relative areas of panel a for a given lambda (l) and range of S  
ep <- function(l, S) {
  
  eps <- 1 - (log(10^(S)))/(4000*l) # Eq. 3
  eps <- ifelse(eps>0,eps,0) # can't have epsilon below 0
  
  areas <- rep(0, length(eps)-1) # initialize vector
  
  if (sum(eps)>0){ # if any feasible values
  
    for (i in 1:length(areas)){
        # subtract area under curve from subsequent area under curve
        # to get area between curves
        areas[i] <- eps[i]-eps[i+1] 
    }

    # relativize to get a probability proxy
    areas <- areas/sum(areas)
  }
  
  return(areas)
}

dfB <- data.frame(matrix(nrow=length(lam), ncol=length(S)-1)) # initialize data frame


for (i in 1:length(lam)) { # run function for array of lambda values
  
  l <- lam[i]
  dfB[i,] <- ep(l, S)
}

colnames(dfB_lines) <- c(paste0('X0', 6:9), paste0('X',10:22))
dfB$lam <- lam

# convert to long
dfB1 <- gather(dfB, key='richness', value='prob', -lam)


# subset the main palette
palB <- pal[2:7]


(pB <- ggplot(dfB1, aes(x=lam, y=prob)) + 
    geom_area(aes(fill=richness), size = 1) +
    scale_fill_manual(values = palB) +
    scale_y_continuous('Probability', 
                       expand = c(0,0)) +
    scale_x_continuous(expression(bold("Speciation rate, "*lambda* " "(My^-1))),
                       # limits=c(0.0035, 0.032),
                       breaks=seq(0.005, 0.03, 0.005),
                       expand = c(0,0)) +
    annotate(geom="text", label='10^6',x=0.011, y = 0.85, color = 'white', size = 5.5, parse=TRUE) +
    annotate(geom="text", label='10^9',x=0.012, y = 0.67, color = 'white', size = 5.5, parse=TRUE) +
    annotate(geom="text", label='10^12',x=0.013, y = 0.52, color = 'white', size = 5.5, parse=TRUE) +
    annotate(geom="text", label='10^15',x=0.014, y = 0.35, color = 'white', size = 5.5, parse=TRUE) +
    annotate(geom="text", label='10^18',x=0.015, y = 0.20, color = 'white', size = 5.5, parse=TRUE) +
    annotate(geom="text", label='10^21',x=0.016, y = 0.07, color = 'white', size = 5.5, parse=TRUE) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 14),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 16, face = "bold"),
          plot.margin = unit(c(1,2,1,1), "lines"))
)

###### panel C #########
# probability of various diversity outcomes across all parameter combinations

probabilities <- function(q, p, MEE=FALSE){
  t <- 4000
  lam <- seq(0.004, 0.030, 0.00001)
  df <- data.frame(lam) # initialize dataframe
  
  pal <- viridis_pal()(7) # create palette
  auc <- c()
  
  index <- 6:23 - 5
  
  for (i in index){
    
    mag <- i + 5
    
    # 10^mag is E(S_t)
    
    if (!MEE) { # if not a mass extinction event, use Eq. 3
      
      rel_ext <- function(lam) 1 - (log(10^mag))/(t*lam)
      
    }else{ # if a mass extinction event, use Eq. 8
      
      rel_ext <- function(lam) 1 - log(10^mag*(1-p)^-1*(1-q*p)^-5)/(lam*3999)
      
    }
    ext <- rel_ext(lam)
    
    # find the x intercept of the function
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
  
  return(area_prob)
}
y0 <- probabilities(0, 1, F) # run function

dfC <- data.frame(mags = 6:22, E_S=y0) 

# plot
pC <- ggplot(dfC, aes(x=mags, y=E_S)) + 
         geom_point(size= 2.0) +
         scale_y_continuous('Probability') +
         scale_x_continuous('Species diversity',
                            breaks=c(6,9, 12, 15, 18, 21),labels = math_format(10^.x)) +
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text = element_text(size = 14),
               axis.line = element_line(colour = "black"),
               axis.title = element_text(size = 16, face = "bold"))


# Save figures and combine
fig1_dir <- paste(figure_dir, "figure1.pdf", sep = "")
plot_grid(pA, pB, pC, align = "hv",labels = "auto", ncol=1)  %>% ggsave2(filename = fig1_dir, units="mm",width = 180, height = 360, dpi = 1000)

