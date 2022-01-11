## halo_tree.R
## Ford Fishman
## 6/30/21

## SETUP ENVIRONMENT
rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)

# packages
library('phytools')
library('seqinr')
library('adephylo')
library('phylosignal')
library('phylobase')
library('ggplot2')
library('ggtree')
library('treeio')
library('dplyr')
library('ggnewscale')
library('RColorBrewer')
library('ggplot2')

read.aln <- read.alignment(file="data/16S/halobacterium_all.afa",format = "fasta")

p.DNAbin <- as.DNAbin(read.aln)

image.DNAbin(p.DNAbin, cex.lab = 0.6)

raxml <- read.raxml("data/16S/tree/halobacterium/RAxML_bipartitionsBranchLabels.T1")
rooted <- root.phylo(raxml@phylo, "NR_113454.1", resolve.root = T)
plot.phylo(rooted, use.edge.length = T,show.node.label = F,show.tip.label = F)

rtt <- distRoot(rooted) # root to tip



# read in ages
ages_file <- paste(data_dir, "16S","halobacterium_ages.csv",sep="/")
ages <- read.csv(ages_file, header = F)
colnames(ages) <- c('acc', 'low_age', 'high_age')
rtt <- rtt[ages$acc]
ages$rtt <- rtt

ggplot(ages, aes(y=rtt, x=low_age)) +
  geom_point() +
  scale_x_continuous('Age (My)') +
  scale_y_log10('Root to Tip Distance') +
  # geom_hline(yintercept = c(5e-4, 1e-4)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

ages2 <- subset(ages, rtt<0.5)


ggplot(ages2, aes(y=rtt, x=low_age)) +
  geom_point() +
  scale_x_continuous('Age (My)') +
  scale_y_log10('Root to Tip Distance') +
  # geom_hline(yintercept = c(5e-4, 1e-4)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"))

summary(lm(rtt~low_age, data=ages2))
