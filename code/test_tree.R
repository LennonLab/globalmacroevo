library('seqinr')
library('ape')

read.aln <- read.alignment(file="./data/p.isolates.afa",format = "fasta")

p.DNAbin <- as.DNAbin(read.aln)

seq.dist.raw <- dist.dna(p.DNAbin, model="raw", pairwise.deletion = F)
nj.tree <- bionj(seq.dist.raw)

outgroup <- match( "Methanosarcina", nj.tree$tip.label)
nj.rooted <- root(nj.tree, outgroup, resolve.root = T)