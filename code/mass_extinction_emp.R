## Modeling impacts of mass extiction on microbial diversity
## Incoroporates estimates of host-associated taxa from EMP
## Ford Fishman

library(ggplot2)
library(tidyr)
library(here)
library(stringr)
library(ggpubr)
library(cowplot)

# Load packages and functions, set environment

setwd(here())
source("code/functions.R")
figure_dir <- "figures/"


# Mass extinction events
GOE <- 2450
OrdovicianSilurian <- 445
LateDevonian <- 375
PermianTriassic <- 252
TriassicJurassic <- 201
CretaceousPaleogene <- 66


EE <- c(OrdovicianSilurian, LateDevonian, PermianTriassic,TriassicJurassic, CretaceousPaleogene)
EE_labels <- c('O-S', 'D', 'P-Tr','Tr-J', 'K-T')
df_EE <- data.frame(year=c(GOE, EE), event=c('GOE',EE_labels))

inEE <- function(mya, EE.=EE){
  return(mya %in% EE.)
}

GOEorEE <- function(mya, EE=EE,GOE=GOE){
  return(mya==GOE|inEE(mya))
}

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
      
      
    } else {
      
      p <- rep(0,3) 
    }
    
    S1[t] <- S1[t-1]*exp(r)*(1-q*p[1])
    S2[t] <- S2[t-1]*exp(r)*(1-q*p[2])
    S3[t] <- S3[t-1]*exp(r)*(1-q*p[3])
    # S4[t] <- S4[t-1]*exp(r)*(1-q*p[4])
  }
  df1 <- data.frame(zero = S1, fifty = S2, ninety = S3, mya=4000:1)
  
  df2 <- gather(df1, key = "intensity", value = "richness",-mya)
  
  df2$intensity <- str_replace(df2$intensity, 'zero', '0%')
  # df2$intensity <- str_replace(df2$intensity, 'ten', '10%')
  df2$intensity <- str_replace(df2$intensity, 'fifty', '50%')
  df2$intensity <- str_replace(df2$intensity, 'ninety', '90%')
  
  return(subset(df2, mya>=10))
}


p1 <- ext_sim(q_min = 1.0) %>%

  ggplot( aes(x = mya, y = log10(richness), group = intensity)) +  
  geom_line(aes(linetype = intensity)) +
  scale_x_continuous("Mya", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 8)) +
  scale_y_continuous("Taxon Diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels = math_format(10^.x), limits = c(0,15)) +
  
  scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed"),
                        labels = c("0%", "50%", "90%")) +
  annotation_logticks(sides = "b") +
  geom_text(data = subset(ext_sim(q_min = 1.0), mya==10),aes(label = intensity), x = Inf, hjust = 1, size =3)+
  geom_text(data=df_EE, aes(x=year,y = 14.5,label = event),size=2, inherit.aes = F)+
  geom_segment(data=df_EE,aes(x = year, xend = year, y = 0, yend=14),size = 4, alpha = 0.4, inherit.aes = F) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),
        plot.margin = unit(c(1,3,1,1), "lines")
  )
p1

p2 <- ext_sim(q_min = 0.1) %>%
  ggplot(aes(x = mya, y = log10(richness), group = intensity)) +  
  geom_line(aes(linetype = intensity)) +
  scale_x_continuous("Mya", trans = reverselog_trans(10), breaks = c(4000, 1000, 100, 10), limits=c(4000, 8)) +
  scale_y_continuous("Taxon Diversity", breaks = c(3, 6, 9, 12), expand = c(0,0),labels = math_format(10^.x), limits = c(0,15)) +
  
  scale_linetype_manual("Taxa removed by each event", 
                        values = c("solid", "dotted", "dashed"),
                        labels = c("0%", "50%", "90%")) +
  annotation_logticks(sides = "b") +
  geom_text(data = subset(ext_sim(q_min = 0.1), mya==10),aes(label = intensity), x = Inf, hjust = 1, size = 3)+
  geom_text(data=df_EE, aes(x=year,y = 14.5,label = event),size=2, inherit.aes = F)+
  geom_segment(data=df_EE,aes(x = year, xend = year, y = 0, yend=14),size = 4, alpha = 0.4, inherit.aes = F) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
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


rel_ext <- function(lam, mag, q, p){
  return(1 - log(10^mag*(1-p)^-1*(1-q*p)^-5)/(lam*3999))
  
}

heatmap_S <- function(q, p){
  lam <- seq(0.002, 0.032, 0.00001)
  df <- data.frame(lam) # initialize dataframe

  pal <- viridis_pal()(7)

  
  for (i in 2:6){
    
    mag <- i*3
    
    # 10^mag is E(S_t)
    ext <- rel_ext(lam, mag, q, p)
    
    # column bind to data frame
    df <- cbind(df,ext)
    names(df)[i] <- paste0("ep",mag)
  }
  
  ep23 <- rel_ext(lam, 23, q, p)
  df <- cbind(df,ep23)
  
  p <- ggplot(df) +
      geom_area(aes(x=lam,y=ep6), fill=pal[2])+
      geom_area(aes(x=lam,y=ep9), fill=pal[3])+
      geom_area(aes(x=lam,y=ep12), fill=pal[4])+
      geom_area(aes(x=lam,y=ep15), fill=pal[5])+
      geom_area(aes(x=lam,y=ep18), fill=pal[6])+
      geom_area(aes(x=lam,y=ep23), fill=pal[7])+
      scale_x_continuous(expression(""*lambda*" (Species/Myr)"),
                         limits=c(0.002, 0.032),
                         breaks=seq(0.005, 0.03, 0.005),
                         expand = c(0,0)) +
      scale_y_continuous("\u03B5", 
                         limits = c(0.0, 0.92),
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
            axis.title = element_text(size = 14, face = "bold"))
  return(p)
}


p3 <- heatmap_S(q=1, p=0.9)
p4 <- heatmap_S(q=0.1, p=0.9)

fig_dir <- paste0(figure_dir, "fig2.png")
plot_grid(p1, p3, p2, p4, align = "hv",labels = "auto", ncol=2, rel_widths = c(1.2,1)) %>% ggsave2(filename = fig_dir, width = 12, height = 8)
