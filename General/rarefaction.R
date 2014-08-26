library(vegan)
library(ggplot2)
library(RColorBrewer)

## Read in total otu counts
otu.counts <- read.table("~/RICE/Lund+Soil/whole/mbn_counts.txt", header = T)
row.names(otu.counts) <- paste("Otu", otu.counts$OTU_ID, sep = "")
otu.counts$OTU_ID <- NULL
map <- read.table("~/RICE/Lund+Soil/soil_lund.map", header = T, row.names = 1)
map <- map[match(colnames(otu.counts), row.names(map)),]
tax <- read.table("~/RICE/Lund+Soil/q.tax", header = T, row.names = 1)
tax <- tax[match(row.names(otu.counts), row.names(tax)),]

## Work on formatting the map
map$BarcodeSequence <- NULL
map$LinkerPrimerSequence <- NULL
map$Field <- NULL
map$Run <- NULL
levels(map$Compartment) <- c(levels(map$Compartment), "Bulk Soil")
map$Compartment[which(map$Compartment == "Bulk_Soil")] <- "Bulk Soil"
map <- subset(map, Compartment != "Pre_planting")
map$Compartment <- factor(map$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
otu.counts <- otu.counts[,match(row.names(map), colnames(otu.counts))]
map$Cultivation <- factor(map$Cultivation, levels = c("EcoFarm", "Organic", "Greenhouse"))
map$Site <- gsub("_", " ", map$Site)
map$Site <- gsub("Sac", "Sacramento", map$Site)
map$Site <- factor(map$Site, levels = c("Arbuckle", "Davis", "Sacramento",
                                        "BB P1", "BB P4", "Ditaler 18", "Ditaler 19",
                                        "DS RR", "Scheidec", "SFT 20 A", "Spooner Airstrip"))
comp.cols <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")

gh.ncounts <- t(subset(t(otu.counts), map$Cultivation == "Greenhouse"))
gh.ncounts <- gh.ncounts[apply(gh.ncounts, 1, sum) > 0,]
gh.map <- subset(map, Cultivation == "Greenhouse")

bs.means <- rowMeans(t(subset(t(gh.ncounts), gh.map$Compartment == "Bulk Soil")))
rs.means <- rowMeans(t(subset(t(gh.ncounts), gh.map$Compartment == "Rhizosphere")))
rp.means <- rowMeans(t(subset(t(gh.ncounts), gh.map$Compartment == "Rhizoplane")))
e.means <- rowMeans(t(subset(t(gh.ncounts), gh.map$Compartment == "Endosphere")))
comp.means <- rbind(bs.means, rs.means, rp.means, e.means)
rc_curve <- function(counts, group, step = 10){
  counts <- round(counts)
  raremax <- min(rowSums(counts))
  created = 0
  for(i in 1:nrow(counts)) {
    j <- 0
    while(j <= raremax){
      value <- rarefy(counts[i,], j)
      sample <- group[i]
      if(created != 1) {
        whole <- c(sample, j, value)
        created = 1
        j = j + step
      } else {
        whole <- rbind(whole, c(sample, j, value))
        j = j + step
      }
    }
  }
  row.names(whole) <- 1:nrow(whole)
  colnames(whole) <- c("Sample", "Samp", "Rich")
  whole <- as.data.frame(whole)
  whole$Rich <- as.numeric(levels(whole$Rich))[whole$Rich]
  whole$Samp <- as.numeric(levels(whole$Samp))[whole$Samp]
  return(whole)
}
gh.curves <- rc_curve(comp.means, c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"), step = 100)
gh.curves$Sample <- factor(gh.curves$Sample, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(gh.curves, aes(x = Samp, y = Rich, color = Sample)) +
  geom_point(alpha = 0.9) +
  theme_bw() +
  labs(x = "Sequences Sampled", y ="Total OTUs Detected", color = "Rhizocompartment") +
  scale_color_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  theme(text = element_text(size = 20), legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape=16, size = 5)))

gh.shan <- cbind(gh.map, Shannon = diversity(t(gh.ncounts)))
ggplot(gh.shan, aes(x = Compartment, y = exp(Shannon), fill = Compartment)) +
  geom_boxplot(alpha = 0.9, outlier.size = 1) +
  scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species", fill = "Rhizocompartment") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20), legend.key = element_blank())




lund.o.rs <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Rhizosphere" & lund.map$Cultivation == "Organic")))
lund.o.rp <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Rhizoplane" & lund.map$Cultivation == "Organic")))
lund.o.e <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Endosphere" & lund.map$Cultivation == "Organic")))
lund.e.rs <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Rhizosphere" & lund.map$Cultivation == "EcoFarm")))
lund.e.rp <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Rhizoplane" & lund.map$Cultivation == "EcoFarm")))
lund.e.e <- rowMeans(t(subset(t(lund.counts), lund.map$Compartment == "Endosphere" & lund.map$Cultivation == "EcoFarm")))
lund.means <- rbind(lund.o.rs, lund.o.rp, lund.o.e, lund.e.rs, lund.e.rp, lund.e.e)
lund.curves <- rc_curve(lund.means, c("Rhizosphere Organic", "Rhizoplane Organic", "Endosphere Organic",
                                      "Rhizosphere EcoFarm", "Rhizoplane EcoFarm", "Endosphere Ecofarm"),
                        step = 100)
lund.curves$Sample <- factor(lund.curves$Sample)


ggplot(lund.curves, aes(x = Samp, y = Rich, color = Sample)) +
  geom_point(alpha = 0.9) +
  theme_bw() +
  labs(x = "Sequences Sampled", y ="Total OTUs Detected", color = "Rhizocompartment") +
  scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#CAB2D6", "#6A3D9A")) + 
  theme(text = element_text(size = 20), legend.key = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape=16, size = 5)))

lund.div <- cbind(lund.map, Shannon = diversity(t(lund.counts)))
lund.div$Compartment <- factor(lund.div$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))
lund.div$Site <- gsub("_", " ", lund.div$Site)
lund.div$Site <- factor(lund.div$Site, levels = c("BB P1", "BB P4", "Ditaler 18", "Ditaler 19",
                                        "DS RR", "Scheidec", "SFT 20 A", "Spooner Airstrip"))
ggplot(lund.div, aes(x = Compartment, y = exp(Shannon), fill = Compartment)) +
  geom_boxplot(alpha = 0.9, outlier.size = 1) +
  scale_fill_manual(values = c("#984EA3", "#4DAF4A", "#377EBA")) +
  facet_grid(.~Site) +
  theme_bw() +
  labs(x = "", y = "Effective Species", fill = "Rhizocompartment") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20))
  
ggplot(lund.div, aes(x = Compartment, y = exp(Shannon), fill = paste(Cultivation, Compartment, sep = " "))) +
  geom_boxplot(alpha = 0.9, outlier.size = 1) +
  scale_fill_manual(values = c("#A6CEE3", "#B2DF8A","#CAB2D6","#1F78B4" ,"#33A02C", "#6A3D9A")) + 
  facet_grid(.~Cultivation) +
  theme_bw() +
  labs(x = "", y = "Effective Species", fill = "Rhizocompartment") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 20))


