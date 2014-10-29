library(ggplot2)
library(reshape2)
library(scales)

setwd("~/RMB/Publication/Data/FieldExp/")

## Load the mapping file
field.map <- read.table("field_map.txt", header = T, row.names = 1)
field.map$BarcodeSequence <- NULL
field.map$LinkerPrimerSequence <- NULL
field.map$Run <- NULL
field.map$Field <- NULL
field.map$lat <- factor(field.map$lat)
field.map$lon <- factor(field.map$lon)

## Load the counts file
field.counts <- read.table("field_counts_with_organelle.txt", header = T, row.names = 1, sep = "\t")
organelle <- read.table("../organelle_otus", header = F)$V1

## Categorize OTUs by microbial or organellar
cats <- rep("Microbial", nrow(field.counts))
cats[row.names(field.counts)%in%organelle] <- "Organellar"

## Sum the counts for the categories
cat.counts <- aggregate(field.counts, as.data.frame(cats), sum)
row.names(cat.counts) <- cat.counts[,1]
cat.counts <- t(cat.counts[,-1])
cat.counts <- cat.counts / rowSums(cat.counts)
cat.counts.m <- melt(cbind(field.map, cat.counts))

## Get summary stats on the counts
cat.sum <- summarySE(cat.counts.m, groupvars = c("variable","Compartment", "Site"), measurevar = "value")
cat.sum$Compartment <- factor(cat.sum$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))


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


