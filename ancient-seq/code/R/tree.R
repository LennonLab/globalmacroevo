## tree.R
## Ford Fishman
## 6/30/21

library('adephylo')
library('ape')
library('phytools')



ml.bootstrap <- read.tree("data/16S/tree/pairs/RAxML_bipartitions.T1")


# outgroup <- match( "MG642124.1", ml.bootstrap$tip.label)
# ml.rooted <- root(ml.bootstrap, outgroup, resolve.root = T)

# ml.rooted$edge.length <- NULL
dists <- distRoot(ml.bootstrap)
print(dists[order(names(dists))])
par(mar=c(1,1,2,1)+0.1)
plotTree(ml.bootstrap,type="fan",fsize=0.7,lwd=1,xlim=c(-1.1,1.1))
# add.scale.bar(cex=0.7)
# nodelabels(ml.bootstrap$node.label, font=2, bg="white", frame="r", cex=0.5)

# x = seq(1,10, 0.1)
# g = function(x) {ifelse(x==10,100,x)}
# y = function(x) {g(x)-3}
# y(x)

