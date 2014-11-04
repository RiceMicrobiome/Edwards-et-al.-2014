library(ggplot2)
library(reshape2)

setwd("~/RMB/Publication/TimeCourse/TC+GH/Data/")


## Load in the files
ts.counts <- read.table("GH_TC_otu_table.txt", header = T, row.names = 1)
ts.map <- read.table("TC+GH.map", header = T, row.names = 1, sep = "\t")
tax <- read.table("TC_GH.tax", header = T, row.names = 1)
mito <- as.character(read.table("mito", header = F)$V1)
plastid <- as.character(read.table("chloro", header = F)$V1)

## Format the data
ts.map.good <- subset(ts.map, Days <= 13)
ts.counts.good <- ts.counts[,match(row.names(ts.map.good), colnames(ts.counts))]
ts.counts.good <- ts.counts.good[,apply(ts.counts.good, 2, sum) > 144]
ts.map.good <- ts.map[match(colnames(ts.counts.good), row.names(ts.map)),]
organelle <- c(mito, plastid)

## Put into DF for platting
ts.sums <- cbind(ts.map.good, RawReads = colSums(ts.counts.good))
ts.sums <- cbind(ts.sums, OrganR = colSums(ts.counts.good[!row.names(ts.counts.good)%in%organelle,]))
ts.sums$Days <- factor(ts.sums$Days)
ts.sums$Rep <- factor(ts.sums$Rep)
ts.sums.m <- melt(ts.sums)
ts.sums.m$variable <- gsub("RawReads", "Raw Reads", ts.sums.m$variable)
ts.sums.m$variable <- gsub("OrganR", "Organellar OTUs Removed", ts.sums.m$variable)
ts.sums.m$Compartment <- factor(ts.sums.m$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
ts.sums.m$variable <- factor(ts.sums.m$variable, levels = c("Raw Reads", "Organellar OTUs Removed"))

## Plot it
ggplot(ts.sums.m, aes(x = Days, y = value, fill = variable)) +
  geom_boxplot() +
  facet_grid(Compartment ~ .) +
  scale_fill_manual(values = c("orange", "darkgreen")) +
  labs(y = "Reads") +
  theme(text = element_text(size = 20))
