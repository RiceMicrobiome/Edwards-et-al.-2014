# Required Packages
library(gplots)
library(ggplot2)
library(plyr)
library(data.table)
library(ape)
library(vegan)
library(VennDiagram)

load("/Users/edwards/RMB/Publication/Data/GreenhouseExp/glm.gh.rda")
tax <- read.table("~/RMB/Publication/Data/GreenhouseExp/gh_tax.txt", header = T, row.names = 1)
thresh.p <- 0.01
da.comp.site.gh.counts <- data.frame(Value = c(
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Arbuckle" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Arbuckle" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Arbuckle" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Arbuckle" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Arbuckle" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Arbuckle" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Davis" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Davis" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Davis" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Davis" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Davis" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Davis" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Sacramento" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Sacramento" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Sacramento" & logFC > 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Endosphere" & Site == "Sacramento" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizoplane" & Site == "Sacramento" & logFC < 0 & padj <= thresh.p)),
  nrow(subset(gh.comp.site.glm, Comp == "Rhizosphere" & Site == "Sacramento" & logFC < 0 & padj <= thresh.p))),
  Site = c(rep("Arbuckle", 6), rep("Davis", 6), rep("Sacramento", 6)),
  Comp= rep(c("Endosphere", "Rhizoplane", "Rhizosphere"), 6),
  Category = rep(c(rep("up" , 3), rep("down", 3)), 3))

da.comp.site.gh.counts$x.vp <- rep(c(rep(7.5 , 3), rep(-8, 3)), 3)
da.comp.site.gh.counts$y.vp <- rep(110, nrow(da.comp.site.gh.counts))
da.comp.site.gh.counts$x.ma <- rep(1, nrow(da.comp.site.gh.counts))
da.comp.site.gh.counts$y.ma <- rep(c(rep(15 , 3), rep(-14, 3)), 3)
da.comp.site.gh.counts$color <- factor(ifelse(da.comp.site.gh.counts$Category == "up", "bs", ifelse(da.comp.site.gh.counts$Comp == "Endosphere", "E", ifelse(da.comp.site.gh.counts$Comp == "Rhizoplane", "RP", "RS"))))

gh.comp.site.glm$sig <- ifelse(gh.comp.site.glm$padj > thresh.p, "ns", "s")
gh.comp.site.glm$color <- ifelse(gh.comp.site.glm$padj > thresh.p, "ns", ifelse(gh.comp.site.glm$logFC < 0, "bs", ifelse(gh.comp.site.glm$Comp == "Endosphere", "E", ifelse(gh.comp.site.glm$Comp == "Rhizoplane", "RP", "RS"))))
gh.comp.site.glm$color <- factor(gh.comp.site.glm$color, levels = c("ns", "bs", "E", "RP", "RS"))

sig.otus <- subset(gh.comp.site.glm, color != "ns")
sig.tax <- tax[tax$OTU%in%unique(sig.otus$OTU),]
sig.otus <- merge(sig.otus, sig.tax, by = "OTU")
write.table(sig.otus, file = "~/RMB/Publication/Data/GreenhouseExp/DA_comp_site_OTUs.txt", sep = "\t")


# Graph MA plots colored on DA OTUs
ggplot(gh.comp.site.glm, aes(x = logCPM, y = logFC, color = color, alpha = sig)) +
  geom_point() +
  scale_color_manual(values = c("grey32", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), guide = FALSE) +
  scale_alpha_manual(values = c(0.2, 1), guide = FALSE) +
  facet_grid(Site ~ Comp) +
  geom_text(aes(x = x.ma, y = y.ma, label = Value, group= NULL), data= da.comp.site.gh.counts, inherit.aes = FALSE, parse = FALSE, size= 8) +
  theme_bw() +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10)) +
  labs(x = "Log 10 Abundance", y = "Log 10 Fold Change")


arb.e.up <- subset(gh.comp.site.glm, color == "E" & Site == "Arbuckle")$OTU
arb.rp.up <- subset(gh.comp.site.glm, color == "RP" & Site == "Arbuckle")$OTU
arb.rs.up <- subset(gh.comp.site.glm, color == "RS" & Site == "Arbuckle")$OTU
dav.e.up <- subset(gh.comp.site.glm, color == "E" & Site == "Davis")$OTU
dav.rp.up <- subset(gh.comp.site.glm, color == "RP" & Site == "Davis")$OTU
dav.rs.up <- subset(gh.comp.site.glm, color == "RS" & Site == "Davis")$OTU
sac.e.up <- subset(gh.comp.site.glm, color == "E" & Site == "Sacramento")$OTU
sac.rp.up <- subset(gh.comp.site.glm, color == "RP" & Site == "Sacramento")$OTU
sac.rs.up <- subset(gh.comp.site.glm, color == "RS" & Site == "Sacramento")$OTU

site.e.up <- list(Arbuckle = arb.e.up, Davis = dav.e.up, Sacramento = sac.e.up)
site.rp.up <- list(Arbuckle = arb.rp.up, Davis = dav.rp.up, Sacramento = sac.rp.up)
site.rs.up <- list(Arbuckle = arb.rs.up, Davis = dav.rs.up, Sacramento = sac.rs.up)
venn.diagram(site.e.up, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/e_up.comp.venn.tiff", cex = 3)
venn.diagram(site.rp.up, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/rp_up.comp.venn.tiff", cex = 3)
venn.diagram(site.rs.up, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/rs_up.comp.venn.tiff", cex = 3)

### Compare statistical significance of compartment overlaps between sites for enriched OTUs
## Endosphere
e_c <- length(unique(c(arb.e.up, dav.e.up, sac.e.up)))
ase <- length(intersect(arb.e.up, sac.e.up))
ade <- length(intersect(arb.e.up, dav.e.up))
sde <- length(intersect(sac.e.up, dav.e.up))
p.ase <- 1 - phyper(ase, length(arb.e.up), e_c - length(arb.e.up), length(sac.e.up))
p.ade <- 1 - phyper(ade, length(arb.e.up), e_c - length(arb.e.up), length(dav.e.up))
p.sde <- 1 - phyper(sde, length(sac.e.up), e_c - length(sac.e.up), length(dav.e.up))

## Rhizoplane
rp_c <- length(unique(c(arb.rp.up, dav.rp.up, sac.rp.up)))
asrp <- length(intersect(arb.rp.up, sac.rp.up))
adrp <- length(intersect(arb.rp.up, dav.rp.up))
sdrp <- length(intersect(sac.rp.up, dav.rp.up))
(p.asrp <- 1 - phyper(asrp, length(arb.rp.up), rp_c - length(arb.rp.up), length(sac.rp.up)))
(p.adrp <- 1 - phyper(adrp, length(arb.rp.up), rp_c - length(arb.rp.up), length(dav.rp.up)))
(p.sdrp <- 1 - phyper(sdrp, length(sac.rp.up), rp_c - length(sac.rp.up), length(dav.rp.up)))

## Rhizosphere
rs_c <- length(unique(c(arb.rs.up, dav.rs.up, sac.rs.up)))
asrs <- length(intersect(arb.rs.up, sac.rs.up))
adrs <- length(intersect(arb.rs.up, dav.rs.up))
sdrs <- length(intersect(sac.rs.up, dav.rs.up))
(p.asrs <- 1 - phyper(asrs, length(arb.rs.up), rs_c - length(arb.rs.up), length(sac.rs.up)))
(p.adrs <- 1 - phyper(adrs, length(arb.rs.up), rs_c - length(arb.rs.up), length(dav.rs.up)))
(p.sdrs <- 1 - phyper(sdrs, length(sac.rs.up), rs_c - length(sac.rs.up), length(dav.rs.up)))

### Depleted OTUs
arb.e.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Arbuckle" & Comp == "Endosphere")$OTU
arb.rp.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Arbuckle" & Comp == "Rhizoplane")$OTU
arb.rs.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Arbuckle" & Comp == "Rhizosphere")$OTU
dav.e.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Davis" & Comp == "Endosphere")$OTU
dav.rp.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Davis" & Comp == "Rhizoplane")$OTU
dav.rs.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Davis" & Comp == "Rhizosphere")$OTU
sac.e.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Sacramento" & Comp == "Endosphere")$OTU
sac.rp.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Sacramento" & Comp == "Rhizoplane")$OTU
sac.rs.down <- subset(gh.comp.site.glm, color == "bs" & Site == "Sacramento" & Comp == "Rhizosphere")$OTU

site.e.down <- list(Arbuckle = arb.e.down, Davis = dav.e.down, Sacramento = sac.e.down)
site.rp.down <- list(Arbuckle = arb.rp.down, Davis = dav.rp.down, Sacramento = sac.rp.down)
site.rs.down <- list(Arbuckle = arb.rs.down, Davis = dav.rs.down, Sacramento = sac.rs.down)
venn.diagram(site.e.down, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/e_down.comp.venn.tiff", cex = 3)
venn.diagram(site.rp.down, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/rp_down.comp.venn.tiff", cex = 3)
venn.diagram(site.rs.down, fill = c("red", "darkorange", "dodgerblue"), filename = "~/RICE/Publication/GHFigures/DA/Site_Comp_Venn/rs_down.comp.venn.tiff", cex = 3)
