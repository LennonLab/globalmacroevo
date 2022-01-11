## SETUP ENVIRONMENT

library('ape')
library('adephylo')
library('phytools')
library('seqinr')
library('treeio')

rm(list=ls()) # removes all objects in the given environment
wd <- "~/GitHub/MicroSpeciation"
data_dir <- paste(wd, "/data/", sep = "")
figure_dir <- paste(wd, "/figures/", sep = "")
getwd()
setwd(wd)

read.aln <- read.alignment(file="data/16S/pairs.afa",format = "fasta")

p.DNAbin <- as.DNAbin(read.aln)

image.DNAbin(p.DNAbin, cex.lab = 0.6)

raxml <- read.raxml("data/16S/tree/pairs/RAxML_bipartitionsBranchLabels.T1")
plotTree(raxml@phylo)

rtt <- distRoot(raxml@phylo)

df <- read.csv("data/16S/pairs.csv")

ancient <- 0
present <- 0

for (i in 1:nrow(df)){
  
  a_n <- df[i,"ancient"]
  a_i <- unname(rtt[a_n])
  p_n <- df[i,"present"]
  p_i <- unname(rtt[p_n])
  
  if (a_i > p_i){
    ancient <- ancient + 1
  } else {
    present <- present + 1
    print(df[i,"ï..taxa"])
  }
}

print(ancient)
