library(ggplot2)
library(RColorBrewer)
library(reshape)
library(VennDiagram)

setwd("~/RMB/Publication/Data/FieldExp/")
load("lund.comp.site.da.rda")

field.tax <- read.table("field_tax.txt", header = T, row.names = 1)

thresh.p <- 0.01
comp.site.da$sig <- ifelse(comp.site.da$padj > thresh.p, "ns", "sig")
comp.site.da$color <- rep("ns", nrow(comp.site.da))
comp.site.da$color[with(comp.site.da, which(sig == "sig" & logFC < 0))] <- "RS"
comp.site.da$color[with(comp.site.da, which(sig == "sig" & Comp == "Endosphere" & logFC > 0))] <- "E"
comp.site.da$color[with(comp.site.da, which(sig == "sig" & Comp == "Rhizoplane" & logFC > 0))] <- "RP"

## Get the factors in the correct order
comp.site.da$Site <- gsub("BBP1", "BB P1", comp.site.da$Site)
comp.site.da$Site <- gsub("BBP4", "BB P4", comp.site.da$Site)
comp.site.da$Site <- gsub("Ditaler18", "Ditaler 18", comp.site.da$Site)
comp.site.da$Site <- gsub("Ditaler19", "Ditaler 19", comp.site.da$Site)
comp.site.da$Site <- gsub("DSRR", "DS RR", comp.site.da$Site)
comp.site.da$Site <- gsub("SFT20A", "SFT 20 A", comp.site.da$Site)
comp.site.da$Site <- factor(comp.site.da$Site, levels = c("BB P1", "BB P4", "Ditaler 18", "Ditaler 19", "Scheidec",
                                                          "DS RR", "SFT 20 A", "Spooner Airstrip"))
comp.site.da$Comp <- factor(comp.site.da$Comp, levels = c("Rhizoplane", "Endosphere"))


## DA Counts
da.counts <- data.frame(Values = c(
  nrow(subset(comp.site.da, color == "E" & Site == "BB P1")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "BB P1")),
  nrow(subset(comp.site.da, color == "RP" & Site == "BB P1")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "BB P1")),
  nrow(subset(comp.site.da, color == "E" & Site == "BB P4")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "BB P4")),
  nrow(subset(comp.site.da, color == "RP" & Site == "BB P4")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "BB P4")),
  nrow(subset(comp.site.da, color == "E" & Site == "Ditaler 18")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "Ditaler 18")),
  nrow(subset(comp.site.da, color == "RP" & Site == "Ditaler 18")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "Ditaler 18")),
  nrow(subset(comp.site.da, color == "E" & Site == "Ditaler 19")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "Ditaler 19")),
  nrow(subset(comp.site.da, color == "RP" & Site == "Ditaler 19")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "Ditaler 19")),
  nrow(subset(comp.site.da, color == "E" & Site == "Scheidec")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "Scheidec")),
  nrow(subset(comp.site.da, color == "RP" & Site == "Scheidec")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "Scheidec")),
  nrow(subset(comp.site.da, color == "E" & Site == "DS RR")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "DS RR")),
  nrow(subset(comp.site.da, color == "RP" & Site == "DS RR")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "DS RR")),
  nrow(subset(comp.site.da, color == "E" & Site == "SFT 20 A")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "SFT 20 A")),
  nrow(subset(comp.site.da, color == "RP" & Site == "SFT 20 A")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "SFT 20 A")),
  nrow(subset(comp.site.da, color == "E" & Site == "Spooner Airstrip")),
  nrow(subset(comp.site.da, Comp == "Endosphere" & color == "RS" & Site == "Spooner Airstrip")),
  nrow(subset(comp.site.da, color == "RP" & Site == "Spooner Airstrip")),
  nrow(subset(comp.site.da, Comp == "Rhizoplane" & color == "RS" & Site == "Spooner Airstrip"))),
  Comp = rep(c("Endosphere", "Endosphere", "Rhizoplane", "Rhizoplane"), 8),
  Site = c(rep("BB P1" , 4), rep("BB P4", 4), rep("Ditaler 18", 4), rep("Ditaler 19", 4), rep("Scheidec", 4), rep("DS RR", 4), rep("SFT 20 A", 4), rep("Spooner Airstrip", 4)))
da.counts$X = rep(-1, nrow(da.counts))
da.counts$Y = c(rep(c(12, -12), 16))

ggplot(comp.site.da, aes(x = logCPM, y = logFC, color = color, alpha = sig)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#377EB8", "grey32", "#4DAF4A", "#984EA3"), guide = FALSE) +
  scale_alpha_manual(values = c(0.2, 1), guide = FALSE) +
  facet_grid(Site ~ Comp) +
  geom_text(aes(x = X, y = Y, label = Values, group= NULL), data= da.counts, inherit.aes = FALSE, parse = FALSE, size= 8) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10)) +
  labs(x = "Log10 Average Abundance", y = "Log10 Fold Change")

## Endosphere enriched OTUs
e.enriched <- subset(comp.site.da, color == "E")
e.enriched.tax <- cbind(e.enriched, field.tax[match(e.enriched$OTU, row.names(field.tax)),])

e.by.site <- as.list(tapply(e.enriched$OTU, e.enriched$Site, unique))
e.all.sites <- Reduce(intersect, e.by.site)
e.all.site.tax <- field.tax[match(e.all.sites, row.names(field.tax)),]

rp.enriched <- subset(comp.site.da, color == "RP")
rp.enriched.tax <- cbind(rp.enriched, field.tax[match(rp.enriched$OTU, row.names(field.tax)),])
rp.by.site <- as.list(tapply(rp.enriched$OTU, rp.enriched$Site, unique))
rp.all.sites <- Reduce(intersect, rp.by.site)

rp.e.enriched.tax <- rbind(rp.enriched.tax, e.enriched.tax)
rp.e.top <- names(sort(table(rp.e.enriched.tax$Phylum), decreasing = T)[1:15])
rp.e.top.tax <- rp.e.enriched.tax[rp.e.enriched.tax$Phylum%in%rp.e.top,]

colors <- c("mediumpurple3", "palegreen4", "skyblue1", "darkorange", "blue", "darkseagreen", "chocolate4",
            "rosybrown1", "firebrick4", "red", "dodgerblue4", "orange", "forestgreen", "grey", "orchid")
ggplot(rp.e.top.tax, aes(x = factor(""), fill = Phylum)) +
  geom_bar(width = 1, position = "fill") +
  coord_polar(theta = "y") +
  theme_bw() +
  scale_fill_manual(values = colors) +
  facet_grid(Comp ~ Site) +
  labs(x = "", y = "") +
  theme(text = element_text(size = 20), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank())

rp.e.enriched.tax$Exp <- rep("Field", nrow(rp.e.enriched.tax))
write.table(rp.e.enriched.tax, file = "~/RICE/Publication/Data/FieldExp/enriched_otus.txt")
rp.e.enriched <- read.csv("~/RMB//Publication//Data//FieldExp//enriched_otus.csv", header = T, row.names = 1)
bbp1.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "BB P1")$OTU)
bbp4.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "BB P4")$OTU)
ditaler18.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "Ditaler 18")$OTU)
ditaler19.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "Ditaler 19")$OTU)
dsrr.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "DS RR")$OTU)
scheidec.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "Scheidec")$OTU)
sfta.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "SFT 20 A")$OTU)
spooner.e <- as.character(subset(rp.e.enriched, color == "E" & Site == "Spooner Airstrip")$OTU)

e.site <- list(bbp1 = bbp1.e, bbp4 = bbp4.e, dt18 = ditaler18.e, 
               dt19 = ditaler19.e, dsr = dsrr.e, sch = scheidec.e,
               sft = sfta.e, spoon = spooner.e)
field.core.e <- field.tax[match(Reduce(intersect, e.site), row.names(field.tax)),]
ggplot(field.core.e, aes(x = Phylum, fill = Class)) +
  geom_bar() +
  coord_flip() +
  theme(text = element_text(size = 30))

bbp1.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "BB P1")$OTU)
bbp4.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "BB P4")$OTU)
ditaler18.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "Ditaler 18")$OTU)
ditaler19.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "Ditaler 19")$OTU)
dsrr.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "DS RR")$OTU)
scheidec.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "Scheidec")$OTU)
sfta.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "SFT 20 A")$OTU)
spooner.rp <- as.character(subset(rp.e.enriched, color == "RP" & Site == "Spooner Airstrip")$OTU)
rp.site <- list(bbp1 = bbp1.rp, bbp4 = bbp4.rp, dt18 = ditaler18.rp, 
               dt19 = ditaler19.rp, dsr = dsrr.rp, sch = scheidec.rp,
               sft = sfta.rp, spoon = spooner.rp)
field.core.rp <- Reduce(intersect, rp.site)
