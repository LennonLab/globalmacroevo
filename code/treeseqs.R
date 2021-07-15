metadata <- read.csv("data/16S/taxonomy.csv", header = F, row.names = 1)
subrates <- read.csv("data/16S/subrates.csv")

phyla <- rep("", nrow(subrates))

for (i in 1:nrow(subrates)){
  acc <- subrates$acc[i]
  phyla[i] <- metadata[acc, 3]
  
}

subrates2 <- cbind(subrates, phyla)

## remove certain rows
fungi <- c("Gemmatimonadota", "Amorphea")

remove <- !( is.na(phyla))
remove <- ifelse(phyla %in% fungi, FALSE, remove)

subrates3 <- subset(subrates2, remove)
subrates4 <- subset(subrates3, low_sub!=0 & high_sub!=0)

write.csv(subrates4, file='data/16S/subrates_pruned.csv',row.names = F)
