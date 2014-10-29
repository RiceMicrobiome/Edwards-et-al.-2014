library(ggplot2)
library(reshape2)
library(scales)

setwd("~/RMB/Publication/Data/GreenhouseExp/")

## Load the mapping file
gh.map <- read.table("gh_map.txt", header = T, row.names = 1)
gh.map$BarcodeSequence <- NULL
gh.map$LinkerPrimerSequence <- NULL
gh.map$Run <- NULL
gh.map$Field <- NULL

## Load the counts file
gh.counts <- read.table("gh_counts_with_organelle.txt", header = T, row.names = 1, sep = "\t")
organelle <- read.table("../organelle_otus", header = F)$V1

## Categorize OTUs by microbial or organellar
cats <- rep("Microbial", nrow(gh.counts))
cats[row.names(gh.counts)%in%organelle] <- "Organellar"

## Sum the counts for the categories
cat.counts <- aggregate(gh.counts, as.data.frame(cats), sum)
row.names(cat.counts) <- cat.counts[,1]
cat.counts <- t(cat.counts[,-1])
cat.counts <- cat.counts / rowSums(cat.counts)
cat.counts.m <- melt(cbind(gh.map, cat.counts))

## Get summary stats on the counts
cat.sum <- summarySE(cat.counts.m, groupvars = c("variable","Compartment", "Site"), measurevar = "value")
cat.sum.cult <- summarySE(cat.counts.m, groupvars = c("variable","Compartment", "Site", "Cultivar"), measurevar = "value")
cat.sum$Compartment <- factor(cat.sum$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
cat.sum.cult$Compartment <- factor(cat.sum.cult$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
cat.sum.cult$Compartment <- factor(cat.sum$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))


## Plot it
ggplot(cat.sum, aes(x = Compartment, y = value, fill = variable, 
                    label = paste(round(value * 100, 1), "%", sep = ""))) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymax = value + se, ymin = value - se), position=position_dodge(0.9), width = 0.2) +
  #geom_text(postion = position_dodge(1), aes(y = 0), vjust = .80) +
  facet_grid(Site ~ .) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  theme_classic() +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "", y = "Percent of Reads", fill = "Read Type") +
  theme(text = element_text(size = 20))

ggplot(cat.sum.cult, aes(x = Compartment, y = value, fill = variable, 
                    label = paste(round(value * 100, 1), "%", sep = ""))) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymax = value + se, ymin = value - se), position=position_dodge(0.9), width = 0.2) +
  #geom_text(postion = position_dodge(1), aes(y = 0), vjust = .80) +
  facet_grid(Site ~ Cultivar) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  theme_classic() +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "", y = "Percent of Reads", fill = "Read Type") +
  theme(text = element_text(size = 20))
