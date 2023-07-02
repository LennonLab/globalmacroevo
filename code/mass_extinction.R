## mass_extinction.R
## Ford Fishman

## Modeling impacts of mass extiction on microbial diversity
## Incoroporates estimates of host-associated taxa from EMP

library(ggplot2)
library(tidyr)
library(here)
library(stringr)
library(cowplot)

# Load packages and functions, set environment
# setwd(here())
source("code/functions.R")
figure_dir <- "figures/"

# Mass extinction events
GOE <- 2450
OrdovicianSilurian <- 445
LateDevonian <- 375
PermianTriassic <- 252
TriassicJurassic <- 201
CretaceousPaleogene <- 66

# create data frame
EE <- c(OrdovicianSilurian, LateDevonian, PermianTriassic,TriassicJurassic, CretaceousPaleogene)
EE_labels <- c('OS', 'D', 'PTr','TrJ', 'KT')
df_EE <- data.frame(year=c(GOE, EE), event=c('GOE',EE_labels))
df_EE$label <- df_EE$year + c(0, 25, -10, 19, -12, 0)

# returns if a timestep is a host-based mass extinction event
inEE <- function(mya, EE.=EE){
  return(mya %in% EE.)
}

# returns if a timestep is a mass extinction event
GOEorEE <- function(mya, EE=EE,GOE=GOE){
  return(mya==GOE|inEE(mya))
}

# simulation function
# q_min = proportion of taxa vulnerable to mass extinction
# lambda = speciation rate
# ep = relative extinction rate (epsilon)
ext_sim <- function(q_min,lambda=0.015,ep=0.5){

  mu <- lambda * ep 
  r <- lambda - mu # diversification rate
  S1 <- c(1) # initial species levels 4Ga - first parameter set
  S2 <- c(1) # initial species levels 4Ga - second parameter set
  S3 <- c(1) # initial species levels 4Ga - third parameter set
  
  for(t in 2:4000) {
    mya <- 4001 - t # how many millions years ago is it?
    q <- 1.0
    
    if (mya==GOE | inEE(mya)) {
      
      p <- c(0, 0.5, 0.9) # proportion of taxa that go extinct
      
      if (mya %in% EE) {
        
        q <- q_min
      }
      
      
    } else { # if not a mass extinction event, set p to 0
      
      p <- rep(0,3) 
    }
    
    S1[t] <- S1[t-1]*exp(r)*(1-q*p[1])
    S2[t] <- S2[t-1]*exp(r)*(1-q*p[2])
    S3[t] <- S3[t-1]*exp(r)*(1-q*p[3])
    # S4[t] <- S4[t-1]*exp(r)*(1-q*p[4])
  }
  # store info in data frame
  df1 <- data.frame(zero = S1, fifty = S2, ninety = S3, mya=4000:1)
  # convert to long format
  df2 <- gather(df1, key = "intensity", value = "richness",-mya)
  
  df2$intensity <- str_replace(df2$intensity, 'zero', '0%')
  df2$intensity <- str_replace(df2$intensity, 'fifty', '50%')
  df2$intensity <- str_replace(df2$intensity, 'ninety', '90%')
  
  return(subset(df2, mya>=10))
}

# annotation styling
params1 <- as.expression(substitute(italic(q) == 1*","~~lambda == 0.015*","~~epsilon == 0.5,))
annotation_p <- as.expression(substitute("Extinction\nintensity ("*italic(p)*")"))

# q = 1
p1 <- ext_sim(q_min = 1.0) %>%

  ggplot( aes(x = mya, y = log10(richness), group = intensity)) +  
  geom_line(aes(linetype = intensity)) +
  scale_x_continuous("Mya", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 8)) +
  scale_y_continuous("Species diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels = math_format(10^.x), limits = c(0,15)) +
  
  scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed"),
                        labels = c("0%", "50%", "90%")) +
  annotation_logticks(sides = "b") +
  geom_text(data = subset(ext_sim(q_min = 1.0), mya==10),aes(label = intensity), x = Inf, hjust = 1, size =4)+
  geom_text(data=df_EE, aes(x=label,y = 14.5,label = event),size=3.5, inherit.aes = F)+
  annotate(geom="text",x=12, y=14.2, label=annotation_p, hjust = 0)+
  geom_segment(data=df_EE,aes(x = year, xend = year, y = 0, yend=14),size = 4, alpha = 0.4, inherit.aes = F) +
  coord_cartesian(xlim = c(4000,8), ylim = c(0, 15), clip = "off") +
  ggtitle(params1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),
        plot.margin = unit(c(1,3,1,1), "lines")
  )
p1

params2 <- as.expression(substitute(italic(q) == 0.1*","~~lambda == 0.015*","~~epsilon == 0.5,))

# q = 0.1

p2 <- ext_sim(q_min = 0.1) %>%
  ggplot(aes(x = mya, y = log10(richness), group = intensity)) +  
  geom_line(aes(linetype = intensity)) +
  scale_x_continuous("Mya", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 8)) +
  scale_y_continuous("Species diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels = math_format(10^.x), limits = c(0,15)) +
  
  scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed"),
                        labels = c("0%", "50%", "90%")) +
  annotation_logticks(sides = "b") +
  geom_text(data = subset(ext_sim(q_min = 0.1), mya==10),aes(label = intensity, y = c(13.1,12.5,11.8)),x=Inf, hjust = 1,size = 4 )+
  geom_text(data=df_EE, aes(x=label,y = 14.5,label = event),size=3.5, inherit.aes = F)+
  annotate(geom="text",x=12, y=14.2, label=annotation_p, hjust = 0)+
  geom_segment(data=df_EE,aes(x = year, xend = year, y = 0, yend=14),size = 4, alpha = 0.4, inherit.aes = F) +
  coord_cartesian(xlim = c(4000,8), ylim = c(0, 15), clip = "off") +
  ggtitle(params2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),
        plot.margin = unit(c(1,3,1,1), "lines")
  )
p2

### heatmaps ###

library('purrr')
library('viridis')


contour_plot <- function(q, p){
  lam <- seq(0.004, 0.030, 0.00001)
  df <- data.frame(lam) # initialize dataframe
  
  pal <- viridis_pal()(8)
  auc <- c()
  
  
  for (i in 2:7){
    
    mag <- i*3
    
    # 10^mag is E(S_t)
    rel_ext <- function(lam) 1 - log(10^mag*(1-p)^-1*(1-q*p)^-5)/(lam*3999)

    ext <- rel_ext(lam)

    
    # column bind to data frame
    df <- cbind(df,ext)
    names(df)[i] <- paste0("ep",mag)
  }
  
  rel_ext <- function(lam) 1 - log(10^23*(1-p)^-1*(1-q*p)^-5)/(lam*3999)
  df <- cbind(df,rel_ext(lam))
  names(df)[ncol(df)] <- paste0("ep","23")
  
  
  params <- as.expression(substitute(italic(p) == a*","~~italic(q) == b,
                                     list(a = format(p),
                                          b = format(q))))
  
  
  
  p <- ggplot(df) +
      geom_area(aes(x=lam,y=ep6), fill=pal[2])+
      geom_area(aes(x=lam,y=ep9), fill=pal[3])+
      geom_area(aes(x=lam,y=ep12), fill=pal[4])+
      geom_area(aes(x=lam,y=ep15), fill=pal[5])+
      geom_area(aes(x=lam,y=ep18), fill=pal[6])+
      geom_area(aes(x=lam,y=ep21), fill=pal[7])+
      geom_area(aes(x=lam,y=ep23), fill=pal[8])+
      scale_x_continuous(expression(bold("Speciation rate, "*lambda* " "(My^-1))),
                         limits=c(0.004, 0.030),
                         breaks=seq(0.005, 0.03, 0.005),
                         expand = c(0,0)) +
      scale_y_continuous("Relative extinction rate, \u03B5", 
                         limits = c(0.0, 1.0),
                         breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), 
                         expand = c(0,0)) +
      ggtitle(params) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill=pal[1]),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 12),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(5,"pt"),
            legend.position = "right",
            legend.text = element_text(size = 12),
            legend.text.align = 1,
            axis.title = element_text(size = 14, face = "bold"))
  return(p)
}

# create plots
p3 <- contour_plot(q=1, p=0.9) + 
  annotate(geom="text", label='10^5',x=0.011, y = 0.58, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^6',x=0.014, y = 0.45, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^9',x=0.0155, y = 0.39, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^12',x=0.0168, y = 0.33, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^15',x=0.018, y = 0.28, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^18',x=0.019, y = 0.23, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^21',x=0.0198, y = 0.19, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^23',x=0.021, y = 0.14, color = 'black', size = 3, parse=TRUE) 
p4 <- contour_plot(q=0.1, p=0.9) +
  annotate(geom="text", label='10^5',x=0.0095, y = 0.71, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^6',x=0.012, y = 0.59, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^9',x=0.0135, y = 0.50, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^12',x=0.015, y = 0.43, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^15',x=0.016, y = 0.37, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^18',x=0.0172, y = 0.31, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^21',x=0.0181, y = 0.26, color = 'white', size = 3, parse=TRUE) +
  annotate(geom="text", label='10^23',x=0.0198, y = 0.19, color = 'black', size = 3, parse=TRUE) 


# combine and save figures
fig_dir <- paste0(figure_dir, "figure2.png")
plot_grid(p1, p3, p2, p4, align = "hv",labels = "auto", ncol=2, rel_widths = c(1.2,1)) %>% ggsave2(filename = fig_dir, width = 12, height = 8.5)
