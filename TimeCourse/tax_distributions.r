library(ggplot2)
library(vegan)
library(reshape2)
library(data.table)
library(scales)

setwd("~/RMB/Publication/TimeCourse/TC+GH/Data/")


#load("GH_TC_norm.rda")
#map <- read.table("TC+GH.map", header = T, row.names = 1, sep = "\t")
#tax <- read.table("TC_GH.tax", header = T, row.names = 1)
#map <- map[match(colnames(tc.gh.norm), row.names(map)),]
#tax <- tax[match(row.names(tc.gh.norm), row.names(tax)),]
#tc.gh.data <- list(tc.gh.norm, map, tax)
#save(tc.gh.data, file = "tc_gh_data.rda")

## Load the data - includes the normalized OTU counts
load("tc_gh_data.rda")
otu.table <- tc.gh.data[[1]]
map <- tc.gh.data[[2]]
tax <- tc.gh.data[[3]]

phy <- aggregate(otu.table, data.frame(tax$Phylum), sum)
row.names(phy) <- phy[,1]
phy <- t(phy[,-1])
phy.top.15.names <- names(sort(colSums(phy), decreasing = T))[1:15]
phy.top.15 <- phy[,colnames(phy)%in%phy.top.15.names]
phy.final <- cbind(phy.top.15, Other = rowSums(phy[,!colnames(phy)%in%phy.top.15.names]))

phy.whole <- cbind(map, phy.final)
phy.whole$Rep <- factor(phy.whole$Rep)
phy.whole$Days <- factor(phy.whole$Days)
phy.whole <- melt(phy.whole)
phy.whole$Compartment <- gsub("Bulk_Soil", "Bulk Soil", phy.whole$Compartment)
phy.whole$Compartment <- factor(phy.whole$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))

phy.col <- c("mediumpurple3", "palegreen4", "skyblue1", "darkorange", "darkseagreen", "chocolate4", "firebrick4", "red", "gold", "goldenrod4", "dodgerblue4", "orange", "forestgreen", "grey", "orchid", "grey50")

ggplot(subset(phy.whole, Cultivation == "Greenhouse"), aes(x = Rep, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Compartment ~ Site) +
  labs(x = "Replicate", y = "") +
  scale_fill_manual(values = phy.col) +
  theme_bw() +
  scale_y_continuous(labels = percent) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20))


ggplot(subset(phy.whole, Cultivation == "TC" & Days != 21 & Days != 34 & Days != 55 & Days != 0), aes(x = Days, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Compartment ~ .) +
  labs(x = "Days", y = "", fill = "Phyla") +
  scale_fill_manual(values = phy.col) +
  theme_bw() +
  scale_y_continuous(labels = percent) +
  theme(text = element_text(size = 20)) +
  guides(fill = guide_legend(reverse = T))


