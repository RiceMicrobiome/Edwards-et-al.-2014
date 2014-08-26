library(ggplot2)
library(RColorBrewer)
library(reshape)
library(VennDiagram)

setwd("~/RICE/Publication/Data/FieldExp/")

load("lund.da.rda")
lund.comp.da <- lund.da[[1]]
lund.cult.da <- lund.da[[2]]
lund.whole <- lund.da[[3]]
field.tax <- read.table("field_tax.txt", header = T, row.names = 1)
lund.map <- lund.whole[,1:7]
lund.map$Field <- NULL
lund.map$Run <- NULL
lund.counts <- t(lund.whole[,8:ncol(lund.whole)])


## Plot differntial MvA plots for differentially abundant OTUs between the compartments
thresh.p <- 0.01
lund.comp.da$sig <- ifelse(lund.comp.da$padj > thresh.p, "ns", "sig")
lund.comp.da$color <- rep("ns", nrow(lund.comp.da))
lund.comp.da$color[with(lund.comp.da, which(sig == "sig" & Comp1 == "Endosphere" & logFC > 0))] <- "E"
lund.comp.da$color[with(lund.comp.da, which(sig == "sig" & Comp1 == "Rhizoplane" & logFC > 0))] <- "RP" 
lund.comp.da$color[with(lund.comp.da, which(sig == "sig" & Comp2 == "Rhizosphere" & logFC < 0))] <- "RS" 
lund.comp.da$color[with(lund.comp.da, which(sig == "sig" & Comp2 == "Rhizoplane" & logFC < 0))] <- "RP"

lund.comp.da.use <- subset(lund.comp.da, Category == "Endosphere-Rhizosphere" | Category == "Rhizoplane-Rhizosphere")
da.counts <- data.frame(Value = c(
  nrow(subset(lund.comp.da.use, color == "E" & Category == "Endosphere-Rhizosphere")),
  nrow(subset(lund.comp.da.use, color == "RS" & Category == "Endosphere-Rhizosphere")),
  nrow(subset(lund.comp.da.use, color == "RP" & Category == "Rhizoplane-Rhizosphere")),
  nrow(subset(lund.comp.da.use, color == "RS" & Category == "Rhizoplane-Rhizosphere"))),
  Category = c("Endosphere-Rhizosphere", "Endosphere-Rhizosphere", "Rhizoplane-Rhizosphere", "Rhizoplane-Rhizosphere"),
  X = c(0, 0, 0, 0),
  Y = c(9, -9, 9, -9))
lund.comp.da.use$Category <- factor(lund.comp.da.use$Category, levels = c("Rhizoplane-Rhizosphere", "Endosphere-Rhizosphere"))
ggplot(lund.comp.da.use, aes(x = logCPM, y = logFC, color = color, alpha = sig)) +
  geom_point() +
  scale_color_manual(values = c("#377EBA", "grey50", "#4DAF4A", "#984EA3"), guide = F) +
  scale_alpha_manual(values = c(0.2, 1), guide = F) +
  facet_grid(.~Category) +
  geom_text(aes(x = X, y = Y, label = Value, group= NULL), data= da.counts, inherit.aes = FALSE, parse = FALSE, size= 8) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x = "Log10 Average Abundance", y = "Log10 Fold Change")


### Get the DA OTUs in rhizoplane and endosphere for making venn diagrams
## Isolate the data
e.up <- subset(lund.comp.da.use, color == "E")
e.down <- subset(lund.comp.da.use, Comp1 == "Endosphere" & color == "RS")
rp.up <- subset(lund.comp.da.use, color == "RP")
rp.down <- subset(lund.comp.da.use, Comp1 == "Rhizoplane" & color == "RS")

## Put into lists for making the venn diagrams
up.venn <- list(Endosphere = e.up$OTU, Rhizoplane = rp.up$OTU)
down.venn <- list(Endosphere = e.down$OTU, Rhizoplane = rp.down$OTU)

## Save the venn diagrams
venn.diagram(up.venn, fill = c("#377EB8", "#4DAF4A"), filename = "~/RICE/Publication/LundFigures/DA/comp.rp.e.up.tiff", cex = 3, cat.cex = 0)
venn.diagram(down.venn, fill = c("#377EB8", "#4DAF4A"), filename = "~/RICE/Publication/LundFigures/DA/comp.rp.e.down.tiff", cex = 3, cat.cex = 0)


## Plot the differentially abundant OTUs between cultivation practices for each compartment
lund.cult.da$color <- gsub("O", "Organic", lund.cult.da$color)
lund.cult.da$color <- gsub("EF", "EcoFarm", lund.cult.da$color)
thresh.p <- 0.01
lund.cult.da$sig <- ifelse(lund.cult.da$padj > thresh.p, "ns", "sig")
lund.cult.da$color <- rep("ns", nrow(lund.cult.da))
lund.cult.da$color[with(lund.cult.da, which(sig == "sig" & Cult1 == "EcoFarm" & logFC > 0))] <- "EcoFarm"
lund.cult.da$color[with(lund.cult.da, which(sig == "sig" & Cult == "EcoFarm" & logFC < 0))] <- "EcoFarm"
lund.cult.da$color[with(lund.cult.da, which(sig == "sig" & Cult1 == "Organic" & logFC > 0))] <- "Organic"
lund.cult.da$color[with(lund.cult.da, which(sig == "sig" & Cult == "Organic" & logFC < 0))] <- "Organic"
cult.da <- data.frame(Values = c(
  nrow(subset(lund.cult.da, color == "EcoFarm" & Comp == "Rhizosphere")),
  nrow(subset(lund.cult.da, color == "Organic" & Comp == "Rhizosphere")),
  nrow(subset(lund.cult.da, color == "EcoFarm" & Comp == "Rhizoplane")),
  nrow(subset(lund.cult.da, color == "Organic" & Comp == "Rhizoplane")),
  nrow(subset(lund.cult.da, color == "EcoFarm" & Comp == "Endosphere")),
  nrow(subset(lund.cult.da, color == "Organic" & Comp == "Endosphere"))),
  Comp = c("Rhizosphere", "Rhizosphere", "Rhizoplane", "Rhizoplane", "Endosphere", "Endosphere"),
  X = c(0, 0, 0, 0, 0, 0),
  Y = c(11, -11, 11, -11, 11, -11))

ggplot(lund.cult.da, aes(x = logCPM, y = logFC, color = color, alpha = sig)) +
  geom_point() +
  scale_color_manual(values = c("brown4", "grey50","darkorange")) +
  scale_alpha_manual(values = c(0.2, 1), guide = F) +
  facet_grid(.~Comp) +
  geom_text(aes(x = X, y = Y, label = Values, group= NULL), data= cult.da, inherit.aes = FALSE, parse = FALSE, size= 8) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.key = element_blank()) +
  labs(x = "Log10 Average Abundance", y = "Log10 Fold Change", color = "Cultivation") +
  guides(colour = guide_legend(override.aes = list(size=5)))

## Join taxonomies with significantly different OTUs
cult.sig <- subset(lund.cult.da, sig == "sig")
cult.sig.tax <- merge(data.table(cult.sig, key = "OTU"), data.table(field.tax, key = "OTU"))




e.e <- subset(lund.cult.da, Comp == "Endosphere" & color == "EcoFarm")
e.e.tax <- cbind(field.tax[match(e.e$OTU, row.names(field.tax)),], e.e)
rp.e <- subset(lund.cult.da, Comp == "Rhizoplane" & color == "EcoFarm")
rp.e.tax <- cbind(field.tax[match(rp.e$OTU, row.names(field.tax)),], rp.e)
rs.e <- subset(lund.cult.da, Comp == "Rhizosphere" & color == "EcoFarm")
rs.e.tax <- cbind(field.tax[match(rs.e$OTU, row.names(field.tax)),], rs.e)

all.eco.tax <- rbind(e.e.tax, rp.e.tax, rs.e.tax)

length(intersect(e.e.top, intersect(rs.e.top, rp.e.top)))
e.otus <- unique(rbind(e.e, rp.e, rs.e)$OTU)
top <- tax[row.names(tax)%in%intersect(e.e.top, intersect(rs.e.top, rp.e.top)),]

e.o <- subset(lund.cult.da, Comp == "Endosphere" & color == "Organic")
e.o.tax <- cbind(field.tax[match(e.o$OTU, row.names(field.tax)),], e.o)
rp.o <- subset(lund.cult.da, Comp == "Rhizoplane" & color == "Organic")
rp.o.tax <- cbind(field.tax[match(rp.o$OTU, row.names(field.tax)),], rp.o)
rs.o <- subset(lund.cult.da, Comp == "Rhizosphere" & color == "Organic")
rs.o.tax <- cbind(field.tax[match(rs.o$OTU, row.names(field.tax)),], rs.o)

all.org.tax <- rbind(e.o.tax, rp.o.tax, rs.o.tax)

comp.cult.tax <- rbind(e.e.tax, rp.e.tax, rs.e.tax, e.o.tax, rp.o.tax, rs.o.tax)
good.phy <- names(sort(table(comp.cult.tax$Phylum), decreasing = T)[1:15])
good.comp.cult.tax <- comp.cult.tax[comp.cult.tax$Phylum%in%good.phy,]
ggplot(good.comp.cult.tax, aes(x = factor(""), fill = Phylum)) +
  geom_bar(stat = "bin", position = "fill", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(Cultivation~Comp) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size = 20), axis.ticks.y = element_blank())

## Proteobacterial Classes
ggplot(subset(good.comp.cult.tax, Class == "Deltaproteobacteria"), aes(x = factor(""), fill = Family)) +
  geom_bar(stat = "bin", position = "fill", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(Cultivation~Comp) +
  #scale_fill_manual(values = colors) +
  #scale_fill_brewer(palette = "Accent") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size = 20), axis.ticks.y = element_blank())


## Anabaena otus
an.counts <- lund.counts[row.names(lund.counts)%in%subset(good.comp.cult.tax, Genus == "Anabaena")$OTU,]
an.df <- cbind(lund.map, Counts = colSums(an.counts))
ggplot(an.df, aes(x = Compartment, y = Counts, fill = Site)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(.~Cultivation) +
  theme_bw() +
  labs(x = "", y = "Counts of Anabaena", fill = "") +
  theme(text = element_text(size = 20))


e.tax <- cbind(field.tax[match(e.otus, row.names(field.tax)),], Cultivation = "EcoFarm")
o.tax <- cbind(field.tax[match(o.otus, row.names(field.tax)),], Cultivation = "Organic")
cultivation.tax <- rbind(e.tax, o.tax)
good <- names(sort(table(cultivation.tax$Phylum), decreasing = T)[1:15])
good.cultivation <- cultivation.tax[cultivation.tax$Phylum%in%good,]
head(good.cultivation)
colors <- c("mediumpurple3", "palegreen4", "skyblue1","darkorange", "blue","darkseagreen", "chocolate4","rosybrown1", "red", "gold", "dodgerblue4", "orange", "forestgreen", "grey", "orchid")
ggplot(good.cultivation, aes(x = factor(""), fill = Phylum)) +
  geom_bar(stat = "bin", position = "fill", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(.~Cultivation) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size = 20))

ggplot(subset(good.cultivation, Phylum == "Actinobacteria"), aes(x = factor(""), fill = Family)) +
  geom_bar(stat = "bin", position = "fill", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(.~Cultivation) +
 # scale_fill_manual(values = colors) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), text = element_text(size = 20))
### Venn diagrams of cultivation overlaps between compartment
ef.otus <- list(Endosphere = e.e$OTU, Rhizoplane = rp.e$OTU, Rhizosphere = rs.e$OTU)
venn.diagram(ef.otus, fill = c("#377EB8", "#4DAF4A", "#984EA3"), filename = "~/RICE/Publication/LundFigures/DA/ef.comp.tiff", cex= 3)
o.otus <- list(Endosphere = e.o$OTU, Rhizoplane = rp.o$OTU, Rhizosphere = rs.o$OTU)
venn.diagram(o.otus, fill = c("#377EB8", "#4DAF4A", "#984EA3"), filename = "~/RICE/Publication/LundFigures/DA/o.comp.tiff", cex= 3)

tax[row.names(tax)%in%"Otu144821",]


## Plot inidividual OTUs
subset(field.tax, row.names(field.tax) == "Otu3533716")
lund.whole$Compartment <- factor(lund.whole$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))
azo <- melt(cbind(lund.map, lund.whole[,colnames(lund.whole)%in%subset(all.org.tax, Genus == "Azospirillum")$OTU]))

ggplot(azo, aes(x = Compartment, y = value, fill = variable)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(Site~Cultivation, scale = "free_x") +
  labs(x = "", y = "Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size = 20))

ggplot(subset(lund.whole, Compartment == "Rhizosphere"), aes(x = Site, y = OtuNew.CleanUp.ReferenceOTU1210826, fill = Site)) +
  geom_boxplot(outlier.size = 1) +
  facet_grid(.~Cultivation, scale = "free_x") +
  labs(x = "", y = "Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), text = element_text(size = 20), legend.key = element_blank())

# Methanobacterium
# Otu77966
# Otu105600
# Otu84255

# Methanocella
# Otu245210
# Otu253768

# 
cor(lund.whole$OtuNew.CleanUp.ReferenceOTU26391, lund.whole$OtuNew.ReferenceOTU1754)
cor(lund.whole$OtuNew.CleanUp.ReferenceOTU26391, lund.whole$OtuNew.CleanUp.ReferenceOTU1210826)
subset(lund.cult.da, OTU == "Otu625556")$color
subset(all.eco.tax, Family == "Syntrophorhabdaceae")
