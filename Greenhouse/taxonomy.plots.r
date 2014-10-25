library(reshape2)
library(ggplot2)
library(vegan)
library(scales)

### Set the working directory
setwd("~/RMB/Publication/Data/GreenhouseExp/")

### Get the counts map and taxonomy table
gh.counts <- read.table("gh_otu_table.txt", header = T, row.names = 1)
gh.map <- read.table("gh_map.txt", header = T, row.names = 1)
gh.tax <- read.table("gh_tax.txt", header = T, row.names = 1)

### Format the map a little bit
gh.map$Field <- NULL
gh.map$BarcodeSequence <- NULL
gh.map$LinkerPrimerSequence <- NULL
gh.map$Run <- NULL
gh.map$Compartment <- factor(gh.map$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
gh.map$Cultivar <- gsub("Nipponebare", "Nipponbare", gh.map$Cultivar)
############################################
## Phylum Plots
############################################
### Aggregate taxa and counts based on phyla
gh.phy <- aggregate(gh.counts, data.frame(gh.tax$Phylum), sum)
row.names(gh.phy) <- gh.phy[,1]
gh.phy <- gh.phy[,-1]
gh.phy.t <- t(gh.phy)
gh.phy.t <- data.frame(gh.phy.t[match(row.names(gh.map), row.names(gh.phy.t)),])
gh.phy.t.rat <- gh.phy.t / rowSums(gh.phy.t)

### Pairwise t.tests on phyla

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


gh.phy.t.davis <- gh.phy.t.rat[match(row.names(subset(gh.map, Site == "Davis")), row.names(gh.phy.t.rat)),]
gh.phy.t.arbuckle <- gh.phy.t.rat[match(row.names(subset(gh.map, Site == "Arbuckle")), row.names(gh.phy.t.rat)),]
gh.phy.t.sac <- gh.phy.t.rat[match(row.names(subset(gh.map, Site == "Sacramento")), row.names(gh.phy.t.rat)),]

phyla.comp.dav <- cbind(phy.pair.wilcox.test(gh.phy.t.davis, subset(gh.map, Site == "Davis")$Compartment), Site = "Davis")
phyla.comp.arb <- cbind(phy.pair.wilcox.test(gh.phy.t.arbuckle, subset(gh.map, Site == "Arbuckle")$Compartment), Site = "Arbuckle")
phyla.comp.sac <- cbind(phy.pair.wilcox.test(gh.phy.t.sac, subset(gh.map, Site == "Sacramento")$Compartment), Site = "Sacramento")

phyla.comp.whole <- rbind(phyla.comp.dav, phyla.comp.arb, phyla.comp.sac)
write.table(phyla.comp.whole, file = "~/RMB/Publication/Data/GreenhouseExp/phyla_wilcox.txt", sep = "\t", quote = F, row.names = F)

ggplot(phyla.comp.whole, aes(x = Taxa, y = -log(padj), color = Site)) +
  geom_point() +
  geom_hline(y = 2.995) +
  facet_grid(Group1 ~ Group2) +
  theme(axis.text.x = element_text(angle = 90))

### Extract the top15 highest represented phyla
top.15 <- names(head(sort(rowSums(gh.phy), decreasing = T), 15))
other <- colSums(gh.phy[!row.names(gh.phy)%in%top.15,])
gh.phy.15 <- rbind(gh.phy[row.names(gh.phy)%in%top.15,], other = other)
gh.phy.15.prop <- gh.phy.15 / colSums(gh.phy.15)
gh.phy.15.long <- melt(cbind(gh.map, t(gh.phy.15.prop)))

gh.phy.15.sum <- summarySE(subset(gh.phy.15.long, Compartment != "Bulk Soil"), groupvars = c("Site", "Compartment", "variable"),
                          measurevar = "value")


ggplot(gh.phy.15.sum, aes(variable, value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2) +
  facet_grid(Site ~ Compartment) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  scale_y_continuous(label = percent, limits = c(0,1)) +
  labs(x = "", y = "Percent of Microbiome") +
  coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size = 20))

### Melt into a long data frame for plotting
gh.phy.whole <- melt(cbind(gh.map, t(gh.phy.15)))
gh.phy.whole$plot_lab <- rep("BS", nrow(gh.phy.whole))
gh.phy.whole$plot_lab[gh.phy.whole$Compartment == "Rhizosphere"] <- "Rs"
gh.phy.whole$plot_lab[gh.phy.whole$Compartment == "Rhizoplane"] <- "Rp"
gh.phy.whole$plot_lab[gh.phy.whole$Compartment == "Endosphere"] <- "E"
gh.phy.whole$plot_lab <- factor(gh.phy.whole$plot_lab, levels = c("BS", "Rs", "Rp", "E"))
### Define the colors
colors <- c("mediumpurple3", "palegreen4", "skyblue1","darkorange", "darkseagreen", "firebrick4", "red", "gold", "goldenrod4", "dodgerblue4", "orange", "forestgreen", "grey", "orchid", "darkmagenta",  "grey50")

### Plot this thing
ggplot(gh.phy.whole, aes(x = plot_lab, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 30), axis.text.y = element_text(size = 20)) +
  guides(fill = guide_legend(reverse = T))

################################################
## See if there are significant differences between 
## certain phyla in the different compartments
################################################



################################################
## Proteobacteria Plots
################################################
### Subset into proteobacteria
prot.tax <- subset(gh.tax, Phylum == "Proteobacteria")
prot.counts <- gh.counts[match(row.names(prot.tax), row.names(gh.counts)),]
gh.prot <- aggregate(prot.counts, data.frame(prot.tax$Class), sum)
row.names(gh.prot) <- gh.prot[,1]
gh.prot <- gh.prot[,-1]
gh.prot.whole <- melt(cbind(gh.map, t(gh.prot)))

prot.cols <- c(colors[1:6], "grey")
ggplot(gh.prot.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = prot.cols) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))

gh.prot.family <- aggregate(prot.counts, data.frame(prot.tax$Order), sum)
row.names(gh.prot.family) <- gh.prot.family[,1]
gh.prot.family <- gh.prot.family[,-1]
top.15 <- names(head(sort(rowSums(gh.prot.family), decreasing = T), 15))
other <- colSums(gh.prot.family[!row.names(gh.prot.family)%in%top.15,])
gh.prot.15 <- rbind(gh.prot.family[row.names(gh.prot.family)%in%top.15,], other = other)
gh.prot.family.whole <- melt(cbind(gh.map, t(gh.prot.15)))
colors.ord <- c("mediumpurple3", "palegreen4", "skyblue1","darkorange", "darkseagreen", "firebrick4", "red", "gold", "goldenrod4", "dodgerblue4", "orange", "forestgreen", "orchid",  "grey","darkmagenta",  "grey50")

ggplot(gh.prot.family.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = colors.ord) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))

comp.cols <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8")
### Plot individual otus
gh.counts.map <- cbind(gh.map, t(gh.counts))
gh.counts.map$Compartment <- factor(gh.counts.map$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(gh.counts.map, aes(x = Site, y = Otu4356171, fill = Compartment)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = comp.cols) +
  theme(text = element_text(size = 20)) +
  labs(x = "")



adiv <- cbind(gh.map, Shannon = diversity(t(gh.counts)))
adiv$Cultivar <- gsub("_", " ", adiv$Cultivar)
adiv$Cultivar <- factor(adiv$Cultivar, levels = c("Glab B", "Glab E", "93-11", "IR50", "M104", "Nipponbare","Soil"))
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")
ggplot(subset(adiv, Compartment == "Rhizosphere"), aes(x = Cultivar, y = exp(Shannon), fill = Cultivar)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = cult.cols) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 30), legend.key = element_blank())

ggplot(subset(adiv, Compartment == "Rhizoplane"), aes(x = Cultivar, y = exp(Shannon), fill = Cultivar)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = cult.cols) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 30), legend.key = element_blank())

ggplot(subset(adiv, Compartment == "Endosphere"), aes(x = Cultivar, y = exp(Shannon), fill = Cultivar)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = cult.cols) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 30), legend.key = element_blank())


ggplot(adiv, aes(x = Cultivar, y = exp(Shannon), fill = Cultivar)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c(cult.cols, "grey")) +
  facet_grid(Compartment~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 30), legend.key = element_blank())


ggplot(adiv, aes(x = Compartment, y = exp(Shannon), fill = Compartment)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 30), legend.key = element_blank())

## T tests to examine pairwise differences in a-divs between all comapartments of all soils
## Do pairwise t-tests on compartments for certain soils
comp_pair_t_test <- function(x) {
  out_df = data.frame(Compartment1 = NA, Compartment2 = NA, p.value = NA)
  comps <- levels(x$Compartment)
  for (i in 1:length(comps)) {
    comp1 <- comps[i]
    for(j in 1:length(comps)) {
      comp2 <- comps[j]
      p.val <- t.test(exp(x[x$Compartment == comp1,]$Shannon), exp(x[x$Compartment == comp2,]$Shannon))$p.value
      out_df <- rbind(out_df, c(comp1, comp2, p.val))
    }
  }
  out_df <- out_df[complete.cases(out_df),]
  final <- remove.dup(out_df)
  final <- final[final$Compartment1 != final$Compartment2,]
  final <- cbind(final, padj = p.adjust(final$p.value, method = "BH"))
  return(final)
}
arb_comp <- comp_pair_t_test(subset(adiv, Site == "Arbuckle"))
dav_comp <- comp_pair_t_test(subset(adiv, Site == "Davis"))
sac_comp <- comp_pair_t_test(subset(adiv, Site == "Sacramento"))

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

comp_site_comparisons <- comp_site_pair_wilcox_test(adiv)
write.table(comp_site_comparisons, file = "comp_site_adiv.txt", quote = F, row.names = F, sep = "\t")
