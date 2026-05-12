#dir
setwd("work/dir/path")

#
library(readr)

#
centroidsFull <- read.csv("path/to/centroids.csv")
View(centroidsFull)

#
lig <- "2_3_BPG"
modern <- "1A3N"
ancestral <- "AncB"
mutant <- "mut_K82N"

# 
modern_centr <- subset(centroidsFull, protein == modern & ligand == lig)
View(modern_centr)

ancestral_centr <- subset(centroidsFull, protein == ancestral & ligand == lig)
View(ancestral_centr)

mutant_centr <- subset(centroidsFull, protein == mutant & ligand == lig)
View(mutant_centr)

#
write.csv(modern_centr, "modern_centr.csv", row.names = FALSE)
write.csv(ancestral_centr, "ancestral_centr.csv", row.names = FALSE)
write.csv(mutant_centr, "mutant_centr.csv", row.names = FALSE)
















