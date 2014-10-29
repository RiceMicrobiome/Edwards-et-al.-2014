library(ggplot2)
library(scales)
library(vegan)
library(reshape2)

setwd("~/RMB/Publication/TimeCourse/TC+GH/Data/")

otu.tab <- read.table("GH_TC_otu_table.txt", header = T, row.names = 1)
map <- read.table("TC+GH.map", header = T, row.names = 1, sep = "\t")
mito <- read.table("mito", header = F)$V1
chloro <- read.table("chloro", header = F)$V1

otu.tab <- otu.tab[,match(row.names(map), colnames(otu.tab))]

categories <- data.frame(OTU = row.names(otu.tab), Category = "Microbial")
categories$Category <- factor(categories$Category, levels = c("Microbial", "Organellar"))
categories$Category[categories$OTU%in%c(as.character(mito), as.character(chloro))] <- "Organellar"

categories$Category <- factor(categories$Category, levels = c("Microbial", "Plastid", "Mitochondria"))
categories$Category[categories$OTU%in%c(mito, chloro)] <- "Organellar"


cat.count <- aggregate(otu.tab, data.frame(categories$Category), sum)
row.names(cat.count) <- cat.count[,1]
cat.count <- t(cat.count[,-1])
cat.rat <- (cat.count / rowSums(cat.count)) * 100
cat.count <- cbind(map, cat.rat[match(row.names(map), row.names(cat.rat)),])
cat.count$Days <- factor(cat.count$Days)
cat.count$Rep <- factor(cat.count$Rep)
cat.count$Compartment <- factor(cat.count$Compartment, levels = c("Bulk_Soil", "Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
gh.ratio <- melt(subset(cat.count, Cultivation == "Greenhouse"))
tc.ratio <- melt(subset(cat.count, Cultivation == "TC"))

ggplot((gh.ratio), aes(x = Rep, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("brown4", "darkorange", "grey50")) +
  facet_grid(Compartment ~ Site) +
  theme_bw() +
  scale_y_continuous(labels = percent) +
  labs(x = "", y = "", fill = "Category") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20))

endo.tc.ratio <- subset(tc.ratio, Compartment == "Endosphere")
endo.gh.ratio <- subset(gh.ratio, Compartment == "Endosphere")
endo.gh.ratio$Days <- as.numeric(as.character(endo.gh.ratio$Days))
endo.gh.ratio$Days[endo.gh.ratio$Days == 30] <- 42
endo.gh.ratio$Days <- factor(endo.gh.ratio$Days)
## find this function in the enriched_analysis.r script
e.tc.summary <- summarySE(endo.tc.ratio, measurevar="value", groupvars=c("variable", "Days"))
e.gh.summary <- summarySE(subset(endo.gh.ratio, Cultivar == "M104"), measurevar="value", groupvars=c("variable", "Days"))

ggplot(subset(e.tc.summary, Days != 21 & Days != 34 & Days != 55), aes(x = Days, y = value / 100, fill = variable, 
                                                                       label = paste(round(value, 2), "%", sep = ""))) +
  geom_bar(stat = "identity", position = "dodge") + #, color = 'black') +
  geom_errorbar(aes(ymin = (value - se) / 100, ymax = (value + se) / 100), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  geom_text(position = position_dodge(1), aes(y = (value + se) / 100), vjust = -.85) +
  theme_classic() +
  scale_y_continuous(labels = percent) +
  labs(x = "Days", y = "Percent of Reads", fill = "Category") +
  theme(text = element_text(size = 30))

ggplot(subset(e.tc.summary, Days != 21 & Days != 34 & Days != 55), aes(x = Days, y = value / 100, fill = variable, 
                                                                       label = paste(round(value, 1), "%", sep = ""))) +
  geom_bar(stat = "identity", position = "dodge") + #, color = 'black') +
  geom_errorbar(aes(ymin = (value - se) / 100, ymax = (value + se) / 100), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  geom_text(position = position_dodge(1), aes(y = 0), vjust = 1) +
  theme_classic() +
  scale_y_continuous(labels = percent) +
  labs(x = "Days", y = "Percent of Reads", fill = "Category") +
  theme(text = element_text(size = 30))

ggplot(e.gh.summary, aes(x = Days, y = value / 100, fill = variable, label = paste(round(value, 1), "%", sep = ""))) +
  geom_bar(stat = "identity", position = "dodge") + #, color = 'black') +
  geom_errorbar(aes(ymin = (value - se) / 100, ymax = (value + se) / 100), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("dodgerblue4", "darkgoldenrod1")) +
  geom_text(position = position_dodge(1), aes(y = 0), vjust = 1) +  theme_classic() +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "Days", y = "Percent of Reads", fill = "Category") +
  theme(text = element_text(size = 30))

ggplot(subset(endo.tc.ratio, Days != 21 & Days != 34 & Days != 55), aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "bin", position = "fill") +
  scale_fill_manual(values = c("brown4", "darkorange", "grey50")) +
  facet_grid(.~Days) +
  theme_bw() +
  scale_y_continuous(labels = percent) +
  labs(x = "", y = "", fill = "Category") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20))
