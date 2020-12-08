length.g <- 1600 # length of gene
cutoff <- 0.03 # % divergence for new OTU
subs <- length.g * cutoff # number of substitutions needed for new OTU

# 16S substitutions per nucleotide per million years (Kuo and Ochman 2009)
low.K16 <- 0.00013
high.K16 <- 0.00091

# number of 16S substitutions per million years
low.sub <- low.K16*length.g
hi.sub <- high.K16*length.g

# Corresponding divergence time (million years):
div.t.ls <- subs / low.sub # divergence time w/ low substitution rate
div.t.hs <- subs / hi.sub # divergence time w/ high substitution rate

# Corresponding speciation rates (million years):
spec.ls <- 1 / div.t.ls
spec.hs <- 1 / div.t.hs

cat("Speciation rates range from", round(spec.ls, 3), "to", round(spec.hs, 3), "species per million years") 
