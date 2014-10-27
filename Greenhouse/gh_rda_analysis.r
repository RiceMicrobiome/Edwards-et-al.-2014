library(vegan)
library(ggplot2)
library(ape)


setwd("~/RMB/Publication/Data/GreenhouseExp/")

gh.map <- read.table("gh_map.txt", header = T, row.names = 1)
gh.map$BarcodeSequence <- NULL
gh.map$LinkerPrimerSequence <- NULL
gh.map$Field <- factor(gh.map$Field)
gh.map$Cultivar <- gsub("Nipponebare", "Nipponbare", gh.map$Cultivar)
gh.map.nosoil <- subset(gh.map, Compartment != "Bulk Soil")

counts <- read.table("gh_otu_table.txt", header = T, row.names = 1)
counts <- counts[,match(row.names(gh.map), colnames(counts))]

gh.wuf <- read.table("../GreenhouseUniFrac/weighted.gh.unifrac", header = T, row.names = 1)
gh.wuf <- gh.wuf[match(row.names(gh.map), row.names(gh.wuf)), match(row.names(gh.map), colnames(gh.wuf))]
gh.wuf.nosoil <- gh.wuf[match(row.names(gh.map.nosoil), row.names(gh.wuf)), match(row.names(gh.map.nosoil), colnames(gh.wuf))]

gh.uuf <- read.table("../GreenhouseUniFrac/unweighted.gh.unifrac", header = T, row.names = 1)
gh.uuf <- gh.uuf[match(row.names(gh.map), row.names(gh.uuf)), match(row.names(gh.map), colnames(gh.uuf))]
gh.uuf.nosoil <- gh.uuf[match(row.names(gh.map.nosoil), row.names(gh.uuf)), match(row.names(gh.map.nosoil), colnames(gh.uuf))]

## Weighted Unifrac distance based RDA
# Everything
wuf.cap.whole <- capscale(as.dist(gh.wuf) ~ Compartment * Site * Cultivar + Condition(Run), data = gh.map, add = T)
anova(wuf.cap.whole, by = "terms")
wuf.cap.whole.axes <- data.frame(cbind(gh.map, scores(wuf.cap.whole)$sites))
wuf.cap.whole.axes$Compartment <- factor(wuf.cap.whole.axes$Compartment, , levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.whole$CCA$eig / sum(wuf.cap.whole$CCA$eig) * 100
comp.col <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.whole.axes, aes(x = CAP1, y = CAP2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (46.0%)", y = "Constrained PCo2 (11.8%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

# Compartment
wuf.cap.comp <- capscale(as.dist(gh.wuf.nosoil) ~ Compartment + Condition(Site + Cultivar + Run), data = gh.map.nosoil, add = T)
anova(wuf.cap.comp)
wuf.cap.comp.axes <- data.frame(cbind(gh.map.nosoil, scores(wuf.cap.comp)$sites))
wuf.cap.comp.axes$Compartment <- factor(wuf.cap.comp.axes$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.comp$CCA$eig / sum(wuf.cap.comp$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.comp.axes, aes(x = CAP1, y = CAP2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1 (89.7%)", y = "Constrained PCo2 (10.25%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

# Site
wuf.cap.site <- capscale(as.dist(gh.wuf.nosoil) ~ Site + Condition(Compartment + Cultivar + Run), data = gh.map.nosoil, add = T)
anova(wuf.cap.site)
wuf.cap.site.axes <- data.frame(cbind(gh.map.nosoil, scores(wuf.cap.site)$sites))
wuf.cap.site.axes$Compartment <- factor(wuf.cap.site.axes$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.site$CCA$eig / sum(wuf.cap.site$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.site.axes, aes(x = CAP1, y = CAP2, shape = Site)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1 (57.7%)", y = "Constrained PCo2 (42.26%)") +
  theme(text = element_text(size = 30))

# Genotype
wuf.cap.geno <- capscale(as.dist(gh.wuf.nosoil) ~ Cultivar + Condition(Site + Compartment + Run), data = gh.map.nosoil, add = T)
anova(wuf.cap.geno)
wuf.cap.geno.axes <- data.frame(cbind(gh.map.nosoil, scores(wuf.cap.geno)$sites))
wuf.cap.geno.axes$Cultivar <- factor(wuf.cap.geno.axes$Cultivar, levels = c("Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))
percent_explained <- wuf.cap.geno$CCA$eig / sum(wuf.cap.geno$CCA$eig) * 100
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")
ggplot(wuf.cap.geno.axes, aes(x = CAP1, y = CAP2, color = Cultivar)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (57.3%)", y = "Constrained PCo2 (16.2%)") +
  scale_color_manual(values = cult.cols) +
  theme(text = element_text(size = 30))


## Unweighted
# Everything
uuf.cap.whole <- capscale(as.dist(gh.uuf) ~ Compartment + Site + Cultivar + Condition(Run), data = gh.map, add = T)
anova(uuf.cap.whole, by = "terms")
uuf.cap.whole.axes <- data.frame(cbind(gh.map, scores(uuf.cap.whole)$sites))
uuf.cap.whole.axes$Compartment <- factor(uuf.cap.whole.axes$Compartment, , levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- uuf.cap.whole$CCA$eig / sum(uuf.cap.whole$CCA$eig) * 100
comp.col <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8")
ggplot(uuf.cap.whole.axes, aes(x = CAP1, y = CAP2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1 (39.7%)", y = "Constrained PCo2 (35.4%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

# Compartment
uuf.cap.comp <- capscale(as.dist(gh.uuf.nosoil) ~ Compartment + Condition(Site + Cultivar + Run), data = gh.map.nosoil, add = T)
anova(uuf.cap.comp)
uuf.cap.comp.axes <- data.frame(cbind(gh.map.nosoil, scores(uuf.cap.comp)$sites))
uuf.cap.comp.axes$Compartment <- factor(uuf.cap.comp.axes$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- uuf.cap.comp$CCA$eig / sum(uuf.cap.comp$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(uuf.cap.comp.axes, aes(x = CAP1, y = CAP2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1 (91.1%)", y = "Constrained PCo2 (8.9%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

# Site
uuf.cap.site <- capscale(as.dist(gh.uuf.nosoil) ~ Site + Condition(Compartment + Cultivar + Run), data = gh.map.nosoil, add = T)
anova(uuf.cap.site)
uuf.cap.site.axes <- data.frame(cbind(gh.map.nosoil, scores(uuf.cap.site)$sites))
uuf.cap.site.axes$Compartment <- factor(uuf.cap.site.axes$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- uuf.cap.site$CCA$eig / sum(uuf.cap.site$CCA$eig) * 100
ggplot(uuf.cap.site.axes, aes(x = CAP1, y = CAP2, shape = Site)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  labs(x = "Constrained PCo1 (69.7%)", y = "Constrained PCo2 (30.2%)") +
  theme(text = element_text(size = 30))

# Genotype
uuf.cap.geno <- capscale(as.dist(gh.uuf.nosoil) ~ Cultivar + Condition(Site + Compartment + Run), data = gh.map.nosoil, add = T)
anova(uuf.cap.geno)
uuf.cap.geno.axes <- data.frame(cbind(gh.map.nosoil, scores(uuf.cap.geno)$sites))
uuf.cap.geno.axes$Cultivar <- factor(uuf.cap.geno.axes$Cultivar, levels = c("Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))
percent_explained <- uuf.cap.geno$CCA$eig / sum(uuf.cap.geno$CCA$eig) * 100
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")
ggplot(uuf.cap.geno.axes, aes(x = CAP1, y = CAP2, color = Cultivar)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (32.8%)", y = "Constrained PCo2 (19.3%)") +
  scale_color_manual(values = cult.cols) +
  theme(text = element_text(size = 30))
