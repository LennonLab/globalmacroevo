library('seqinr')
library('ape')
library('phylobase')

# read.aln <- read.alignment(file="./data/16S/trios_icecore/test_trio.afa",format = "fasta")
read.aln <- read.alignment(file="./data/16S/trios_icecore/test_trio2.afa",format = "fasta")

p.DNAbin <- as.DNAbin(read.aln)

image.DNAbin(p.DNAbin, cex.lab = 0.6)

(seq.dist.raw <- dist.dna(p.DNAbin, model="raw", pairwise.deletion = F))
nj.tree <- bionj(seq.dist.raw)

# outgroup <- match( "NR_044774.1", nj.tree$tip.label)
outgroup <- match( "NR_121697.2", nj.tree$tip.label)
nj.rooted <- root(nj.tree, outgroup, resolve.root = T)

rooted <- as(nj.rooted, 'phylo4d')

plot.phylo(nj.rooted, use.edge.length = T,show.node.label = T, node.depth = 1)
