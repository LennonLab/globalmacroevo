## probability_curves.R
## Ford Fishman

library(ggplot2)
library(viridis)
library(tidyr)
library(here)

wd <- here()
data_dir <- paste0(wd, "/data/")
figure_dir <- paste0(wd, "/figures/")

lam <- seq(0.002, 0.032, 0.00001)
S <- c(6,9,12,15,18,21,23)

ep <- function(l, S) {
  
  eps <- 1 - (log(10^(S)))/(4000*l)
  eps <- ifelse(eps>0,eps,0)
  
  areas <- rep(0, length(eps)-1)
  
  if (sum(eps)>0){
    
    for (i in 1:length(areas)){
      areas[i] <- eps[i]-eps[i+1]
    }
    areas <- areas/sum(areas)
  }
  
  # return(areas)
  return(areas)
}

df <- data.frame(matrix(nrow=length(lam), ncol=length(S)-1))

# colnames(df) <- S

for (i in 1:length(lam)) {
  
  l <- lam[i]
  df[i,] <- ep(l, S)
}

df$lam <- lam

df1 <- gather(df, key='richness', value='prob', -lam)
annotation <- data.frame(
  x = c(0.0075,4.5),
  y = c(20,25),
  label = c("label 1", "label 2")
)
pal6 <- viridis_pal()(15)[14]
pal <- viridis_pal()(7)[2:7]
pal[6] <- pal6

(p1 <- ggplot(df1, aes(x=lam, y=prob, fill=richness)) + 
  # geom_line(aes(color=richness), size = 2.0, alpha=0.4) +
  geom_area(size = 1) +
  scale_fill_manual(values = pal) +
  scale_y_continuous('Probability', 
                     expand = c(0,0)) +
  scale_x_continuous(expression(bold("Speciation rate, "*lambda*" (Species/Myr)")),
                     limits=c(0.0035, 0.032),
                     breaks=seq(0.005, 0.03, 0.005),
                     expand = c(0,0)) +
  annotate(geom="text", label='10^6',x=0.011, y = 0.85, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^9',x=0.012, y = 0.67, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^12',x=0.013, y = 0.52, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^15',x=0.014, y = 0.35, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^18',x=0.015, y = 0.20, color = 'white', size = 6, parse=TRUE) +
  annotate(geom="text", label='10^21',x=0.016, y = 0.07, color = 'black', size = 6, parse=TRUE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, face = "bold")))

fig2_dir <- paste(figure_dir, "figure2.png", sep = "")
ggsave(plot = p1, filename = fig2_dir, width = 7, height = 5)
