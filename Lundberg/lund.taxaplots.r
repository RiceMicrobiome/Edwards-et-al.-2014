library(reshape2)
library(ggplot2)
library(scales)

### Set the working directory
setwd("~/RMB/Publication/Data/FieldExp/")

### Load in the data
field.counts <- read.table("field_otu_table.txt", header = T, row.names = 1)
field.map <- read.table("field_map.txt", header = T, row.names = 1)
field.tax <- read.table("field_tax.txt", header = T, row.names = 1)

### Format the map a bit
field.map$Field <- NULL
field.map$Run <- NULL
field.map$BarcodeSequence <- NULL
field.map$LinkerPrimerSequence <- NULL

field.counts <- field.counts[,match(row.names(field.map), colnames(field.counts))]

############################################
## Phylum Plots
############################################
### Aggregate taxa and counts based on phyla
field.phy <- aggregate(field.counts, data.frame(field.tax$Phylum), sum)
row.names(field.phy) <- field.phy[,1]
field.phy <- field.phy[,-1]
field.phy.rat <- field.phy / colSums(field.phy)

### Extract the top15 highest represented phyla
top.15 <- names(head(sort(rowSums(field.phy), decreasing = T), 15))
other <- colSums(field.phy[!row.names(field.phy)%in%top.15,])
field.phy.15 <- rbind(field.phy[row.names(field.phy)%in%top.15,], other = other)
field.phy.15.ratio <- field.phy.15 / colSums(field.phy.15)
field.phy.15.ratio <- field.phy.15.ratio[,match(row.names(field.map), names(field.phy.15.ratio))]

### Melt into a long data frame for plotting
field.phy.whole <- melt(cbind((field.map[,-2:-3])[,-4], t(field.phy.15.ratio)))
field.phy.15.sum <- summarySE(field.phy.whole, groupvars = c("Site", "Compartment", "variable"),
                           measurevar = "value")

field.phy.15.sum$Compartment <- factor(field.phy.15.sum$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))


### Define the colors
colors <- c("mediumpurple3", "palegreen4", "skyblue1","darkorange", "blue","darkseagreen", "chocolate4","rosybrown1", "red", "gold", "dodgerblue4", "orange", "forestgreen", "grey", "orchid",  "grey50")

### Plot this thing
ggplot(field.phy.15.sum, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2, position = position_dodge(0.9)) +
  facet_grid(Compartment ~ Site) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_y_continuous(label = percent, limits = c(0,1)) +
  labs(x = "", y = "Percent of Microbiome") +
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size = 20))


ggplot(field.phy.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 30), axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 40, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))

### Get pvalues for significant phyla abundance between compartment
phy.pair.wilcox.test <- function(phyla_counts, groups, method = "BH") {
  phyla_counts <- data.frame(phyla_counts)
  taxa <- names(phyla_counts)
  group_names <- as.character(unique(groups))
  out.df <- data.frame(Taxa = NA, Group1 = NA, Group1_Mean = NA, Group2 = NA, Group2_Mean = NA, pval = NA, padj = NA)
  
  # Loop through taxa
  for (i in 1:ncol(phyla_counts)) {
    tax <- taxa[i]
    tax.counts <- phyla_counts[,i]
    temp.df <- data.frame(Taxa = NA, Group1 = NA, Group1_Mean = NA, Group2 = NA, Group2_Mean = NA, pval = NA)
    for (j in 1:length(group_names)) {
      group1 <- group_names[j]
      group1.counts <- as.numeric(tax.counts[groups == group1])
      group1.mean <- mean(group1.counts)
      #print(c(length(group1.counts), group1))
      for (k in 1:length(group_names)) {
        group2 <- group_names[k]
        group2.counts <- as.numeric(tax.counts[groups == group2])
        group2.mean <- mean(group2.counts)
        #print(c(length(group1.counts), length(group2.counts)))
        p.val <- wilcox.test(group1.counts, group2.counts)$p.value
        temp.df <- rbind(temp.df, c(tax, group1, group1.mean, group2, group2.mean, p.val))
      }
    }
    temp.df <- remove.dup(temp.df)
    temp.df <- temp.df[complete.cases(temp.df),]
    temp.df <- temp.df[temp.df$Group1 != temp.df$Group2,]
    #print(temp.df)
    temp.df$padj <- p.adjust(temp.df$pval, method = method)
    out.df <- rbind(out.df, temp.df)
  }
  out.df <- out.df[complete.cases(out.df),]
  out.df <- remove.dup(out.df)
  return(out.df)
}
field.phy2 <- t(field.phy.rat[,match(row.names(field.map), names(field.phy.rat))])
bbp1.phy <- field.phy2[match(row.names(subset(field.map, Site == "BB P1")), row.names(field.phy2)),]
bbp4.phy <- field.phy2[match(row.names(subset(field.map, Site == "BB P4")), row.names(field.phy2)),]
d18.phy <- field.phy2[match(row.names(subset(field.map, Site == "Ditaler 18")), row.names(field.phy2)),]
d19.phy <- field.phy2[match(row.names(subset(field.map, Site == "Ditaler 19")), row.names(field.phy2)),]
dsrr.phy <- field.phy2[match(row.names(subset(field.map, Site == "DS RR")), row.names(field.phy2)),]
scheidec.phy <- field.phy2[match(row.names(subset(field.map, Site == "Scheidec")), row.names(field.phy2)),]
sft.phy <- field.phy2[match(row.names(subset(field.map, Site == "SFT 20 A")), row.names(field.phy2)),]
spoon.phy <- field.phy2[match(row.names(subset(field.map, Site == "Spooner Airstrip")), row.names(field.phy2)),]

bb1.test <- cbind(phy.pair.wilcox.test(phyla_counts = bbp1.phy, groups = subset(field.map, Site == "BB P1")$Compartment), Site = "BB P1")
bb4.test <- cbind(phy.pair.wilcox.test(phyla_counts = bbp4.phy, groups = subset(field.map, Site == "BB P4")$Compartment), Site = "BB P4")
d18.test <- cbind(phy.pair.wilcox.test(phyla_counts = d18.phy, groups = subset(field.map, Site == "Ditaler 18")$Compartment), Site = "Ditaler 18")
d19.test <- cbind(phy.pair.wilcox.test(phyla_counts = d19.phy, groups = subset(field.map, Site == "Ditaler 19")$Compartment), Site = "Ditaler 19")
dsrr.test <- cbind(phy.pair.wilcox.test(phyla_counts = dsrr.phy, groups = subset(field.map, Site == "DS RR")$Compartment), Site = "DS RR")
scheidec.test <- cbind(phy.pair.wilcox.test(phyla_counts = scheidec.phy, groups = subset(field.map, Site == "Scheidec")$Compartment), Site = "Scheidec")
sft.test <- cbind(phy.pair.wilcox.test(phyla_counts = sft.phy, groups = subset(field.map, Site == "SFT 20 A")$Compartment), Site = "SFT 20 A")
spoon.test <- cbind(phy.pair.wilcox.test(phyla_counts = spoon.phy, groups = subset(field.map, Site == "Spooner Airstrip")$Compartment), Site = "Spooner Airstrip")


field.phy.w.tests <- rbind(bb1.test, bb4.test, d18.test, d19.test,
                           dsrr.test, scheidec.test, sft.test, spoon.test)
field.phy.w.tests$Group1_Mean <- as.numeric(as.character(field.phy.w.tests$Group1_Mean)) * 100
field.phy.w.tests$Group2_Mean <- as.numeric(as.character(field.phy.w.tests$Group2_Mean)) * 100

write.table(field.phy.w.tests, file = "comp_site_phyla_wilcox.txt", row.names = F, sep = "\t", quote = F) 

### Plot Proteobacterial classes
prot.tax <- subset(field.tax, Phylum == "Proteobacteria")
prot.counts <- field.counts[match(row.names(prot.tax), row.names(field.counts)),]
field.prot <- aggregate(prot.counts, data.frame(prot.tax$Class), sum)
row.names(field.prot) <- field.prot[,1]
field.prot <- field.prot[,-1]
field.prot.whole <- melt(cbind(field.map, t(field.prot)))

prot.cols <- c(colors[1:4], colors[6], "grey")
field.prot.whole$Compartment <- factor(field.prot.whole$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(field.prot.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = prot.cols) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))

### Alpha Diversity
comp_site_pair_wilcox_test <- function(x) {
  out_df = data.frame(Compartment1 = NA, Site1 = NA, Compartment2 = NA, Site2 = NA, Group1_Mean = NA, Group2_Mean = NA, p.value = NA)
  x$compsite <- factor(paste(x$Compartment, x$Site, sep = "_"))
  compsites <- levels(x$compsite)
  for (i in 1:length(compsites)) {
    compsite1 <- compsites[i]
    comp1 <- strsplit(compsite1, split = "_")[[1]][1]
    site1 <- strsplit(compsite1, split = "_")[[1]][2]
    mean1 <- mean(exp(x[x$Compartment == comp1 & x$Site == site1,]$Shannon))
    for(j in 1:length(compsites)) {
      compsite2 <- compsites[j]
      comp2 <- strsplit(compsite2, split = "_")[[1]][1]
      site2 <- strsplit(compsite2, split = "_")[[1]][2]
      mean2 <- mean(exp(x[x$Compartment == comp2 & x$Site == site2,]$Shannon))
      p.val <- wilcox.test(exp(x[x$Compartment == comp1 & x$Site == site1,]$Shannon), exp(x[x$Compartment == comp2 & x$Site == site2,]$Shannon))$p.value
      out_df <- rbind(out_df, c(comp1, site1, comp2, site2, mean1, mean2, p.val))
    }
  }
  out_df <- out_df[complete.cases(out_df),]
  final <- remove.dup(out_df)
  final <- final[paste(final$Compartment1, final$Site1) != paste(final$Compartment2, final$Site2),]
  final <- cbind(final, padj = p.adjust(final$p.value, method = "BH"))
  return(final)
}

field.adiv <- cbind(field.map, Shannon = diversity(t(field.counts)))
comp_site_comparisons <- comp_site_pair_wilcox_test(field.adiv)
write.table(comp_site_comparisons, file = "comp_site_adiv_wilcox.txt", sep = "\t", quote = F, row.names = F)
