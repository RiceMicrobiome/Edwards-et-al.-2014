library(ggplot2)
library(reshape2)
library(scales)

setwd("~/RMB/Publication/Data/Dip/")
source("~/RMB/Publication/Scripts/General/rmb_functions.r")

counts <- read.table("dip_otu_table.txt", header = T, row.names = 1)
mito <- read.table("~/RMB/Reference/mito_otus.txt", header = F)$V1
plastid <- read.table("~/RMB/Reference/plastid_otus.txt", header = F)$V1

map <- read.table("dip.map", header = T, row.names = 1, sep = "\t")
map <- map[match(colnames(counts), row.names(map)),]

category <- rep("Microbial", nrow(counts))
category[row.names(counts)%in%c(mito, plastid)] <- "Organellar"


cat.counts <- aggregate(counts, data.frame(category), sum)
row.names(cat.counts) <- cat.counts[,1]
cat.counts <- t(cat.counts[,-1])
cat.ratio <- (cat.counts / rowSums(cat.counts)) * 100
cat.map <- map[match(row.names(cat.ratio), row.names(map)),]
cat.m <- melt(cbind(cat.map, cat.ratio))
cat.m.sum <- summarySE(cat.m, measurevar = "value", groupvars = c("SampleType", "variable"))
cat.m.sum$SampleType <- factor(cat.m.sum$SampleType, levels = c("0h Pre", "0h Post", "24h Post", "Soil"))

ggplot(subset(cat.m.sum, SampleType != "Soil"), aes(x = SampleType, y = value / 100, 
                                                    fill = variable, label = paste(round(value, 2), "%", sep = ""))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (value / 100) - (se / 100), ymax = (value / 100) + (se / 100)), width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_manual(values = c("royalblue4", "orange")) +
  geom_text(position = position_dodge(1), vjust = -0.8) +
  theme_classic() +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "", y = "Percent of Reads", fill = "Read Type") +
  theme(text = element_text(size = 20))

t.test(subset(cat.m, variable == "Microbial" & SampleType == "0h Pre")$value, 
       subset(cat.m, variable == "Microbial" & SampleType == "0h Post")$value)

