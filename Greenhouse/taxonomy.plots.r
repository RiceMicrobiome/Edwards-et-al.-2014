library(reshape2)
library(ggplot2)
library(vegan)

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

############################################
## Phylum Plots
############################################
### Aggregate taxa and counts based on phyla
gh.phy <- aggregate(gh.counts, data.frame(gh.tax$Phylum), sum)
row.names(gh.phy) <- gh.phy[,1]
gh.phy <- gh.phy[,-1]

### Extract the top15 highest represented phyla
top.15 <- names(head(sort(rowSums(gh.phy), decreasing = T), 15))
other <- colSums(gh.phy[!row.names(gh.phy)%in%top.15,])
gh.phy.15 <- rbind(gh.phy[row.names(gh.phy)%in%top.15,], other = other)


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
adiv$Cultivar <- gsub("Nipponebare", "Nipponbare", adiv$Cultivar)
adiv$Cultivar <- factor(adiv$Cultivar, levels = c("Glab B", "Glab E", "93-11", "IR50", "M104", "Nipponbare","Soil"))
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")
ggplot(subset(adiv, Compartment == "Rhizosphere"), aes(x = Cultivar, y = exp(Shannon), fill = Cultivar)) +
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
remove.dup <- function(df) {
  df.sort <- t(apply(df, 1, sort))
  df.sort <- df[!duplicated(df.sort),]
  return(df.sort)
}

comp_pair_t_test <- function(x) {
  out_df = data.frame(Compartment1 = NA, Compartment2 = NA, p.value = NA)
  comps <- levels(x$Compartment)
  for (i in 1:length(comps)) {
    comp1 <- comps[i]
    for(j in 1:length(comps)) {
      comp2 <- comps[j]
      p.val <- t.test(x[x$Compartment == comp1,]$Shannon, x[x$Compartment == comp2,]$Shannon)$p.value
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
