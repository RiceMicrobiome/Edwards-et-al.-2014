library(vegan)
library(ggplot2)
library(ape)
library(RColorBrewer)

setwd("~/RMB/Publication/Data/FieldExp/")

field.map <- read.table("field_map.txt", header = T, row.names = 1)
field.map$BarcodeSequence <- NULL
field.map$LinkerPrimerSequence <- NULL
field.map$Field <- NULL

field.wuf <- read.table("../FieldUniFrac/weighted.field.unifrac", header = T, row.names = 1)
field.wuf <- field.wuf[match(row.names(field.map), row.names(field.wuf)), match(row.names(field.map), colnames(field.wuf))]

field.uuf <- read.table("../FieldUniFrac/unweighted.field.unifrac", header = T, row.names = 1)
field.wuf <- field.uuf[match(row.names(field.map), row.names(field.uuf)), match(row.names(field.map), colnames(field.uuf))]

## Examine Compartment Separation
# Weighted
wuf.cap.comp <- capscale(as.dist(field.wuf) ~ Compartment + Condition(Site + Cultivation + Run), data = field.map, add = T)
anova(wuf.cap.comp)
var_get(wuf.cap.comp)
wuf.cap.comp.axes <- data.frame(cbind(field.map, scores(wuf.cap.comp)$sites))
wuf.cap.comp.axes$Compartment <- factor(wuf.cap.comp.axes$Compartment, , levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.comp$CCA$eig / sum(wuf.cap.comp$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.comp.axes, aes(x = CAP1, y = CAP2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (87.6%)", y = "Constrained PCo2 (12.4%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

# Unweighted
uuf.cap.comp <- capscale(as.dist(field.uuf) ~ Compartment + Condition(Site + Cultivation + Run), data = field.map, add = T)
anova(uuf.cap.comp)
var_get(uuf.cap.comp)
uuf.cap.comp.axes <- data.frame(cbind(field.map, scores(uuf.cap.comp)$sites))
uuf.cap.comp.axes$Compartment <- factor(uuf.cap.comp.axes$Compartment, , levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- uuf.cap.comp$CCA$eig / sum(uuf.cap.comp$CCA$eig) * 100
comp.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(uuf.cap.comp.axes, aes(x = CAP1, y = CAP2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "Constrained PCo1 (83.2%)", y = "Constrained PCo2 (16.8%)") +
  scale_color_manual(values = comp.col) +
  theme(text = element_text(size = 30))

## Site
# Weighted
wuf.cap.site <- capscale(as.dist(field.wuf) ~ Site + Condition(Compartment + Cultivation + Run), data = field.map, add = T)
anova(wuf.cap.site)
var_get(wuf.cap.site)
wuf.cap.site.axes <- data.frame(cbind(field.map, scores(wuf.cap.site)$sites))
wuf.cap.site.axes$siteartment <- factor(wuf.cap.site.axes$Site, , levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
percent_explained <- wuf.cap.site$CCA$eig / sum(wuf.cap.site$CCA$eig) * 100
site.col <- c("#984EA3", "#4DAF4A", "#377EB8")
ggplot(wuf.cap.site.axes, aes(x = CAP1, y = CAP2, color = factor(lat))) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  scale_color_manual(values = brewer.pal(n = 10, "YlGnBu")[3:10]) +
  labs(x = "Constrained PCo1 (53.6%)", y = "Constrained PCo2 (19.5%)", color = "Latitude") +
  theme(text = element_text(size = 30))

# Unweighted
uuf.cap.site <- capscale(as.dist(field.uuf) ~ Site + Condition(Compartment + Cultivation + Run), data = field.map, add = T)
anova(uuf.cap.site)
var_get(uuf.cap.site)
uuf.cap.site.axes <- data.frame(cbind(field.map, scores(uuf.cap.site)$sites))
percent_explained <- uuf.cap.site$CCA$eig / sum(uuf.cap.site$CCA$eig) * 100

ggplot(uuf.cap.site.axes, aes(x = CAP1, y = CAP2, color = factor(lat))) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.75) +
  theme_classic() +
  scale_color_manual(values = brewer.pal(n = 10, "YlGnBu")[3:10]) +
  labs(x = "Constrained PCo1 (47.4%)", y = "Constrained PCo2 (19.8%)", color = "Latitude") +
  theme(text = element_text(size = 30))

#######!!!!!!!!!!!!
# This doesn't work
## Cultivation
# Weighted
wuf.cap.cult <- capscale(as.dist(field.wuf) ~ Cultivation + Condition(Compartment + Site + Run), data = field.map, add = T)
anova(wuf.cap.cult)
wuf.cap.cult.axes <- data.frame(cbind(field.map, scores(wuf.cap.cult)$sites))
percent_explained <- wuf.cap.cult$CCA$eig / sum(wuf.cap.cult$CCA$eig) * 100
ggplot(wuf.cap.cult.axes, aes(x = MDS1, y = MDS2, color = Cultivation)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  scale_color_manual(values = c("brown4", "darkorange")) +
  labs(x = "Constrained PCo1 (53.6%)", y = "Constrained PCo2 (19.5%)") +
  theme(text = element_text(size = 30))

# Unweighted
uuf.cap.cult <- capscale(as.dist(field.uuf) ~ Cultivation + Condition(Compartment + Site + Run), data = field.map, add = T)
anova(uuf.cap.cult)
uuf.cap.cult.axes <- data.frame(cbind(field.map, uuf.cap.site$CCA$wa))
percent_explained <- uuf.cap.cult$CCA$eig / sum(uuf.cap.cult$CCA$eig) * 100
ggplot(uuf.cap.cult.axes, aes(x = CAP1, y = CAP2, color = Cultivation)) +
  geom_vline(x = 0, alpha = 0.9) +
  geom_hline(y = 0, alpha = 0.9) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  scale_color_manual(values = c("brown4", "darkorange")) +
  labs(x = "Constrained PCo1 (53.6%)", y = "Constrained PCo2 (19.5%)") +
  theme(text = element_text(size = 30))


uuf.cap.cult <- capscale(as.dist(field.uuf) ~ Compartment + Cultivation + Condition(Run), data = field.map, add = T)
anova(uuf.cap.cult, by = "terms")
