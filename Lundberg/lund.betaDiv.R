library(vegan)
library(reshape2)
library(ggplot2)
library(ape)


### Set working directory for this analysis
setwd("~/RMB/Publication/Data/FieldUniFrac/")

### Load in the data
# Load in the mapping file and format it a bit
field.map <- read.table("~/RMB/Publication/Data/FieldExp/field_map.txt", header = T, row.names = 1)

# Load in weighted and unweighted unifrac distances
wuf.field.df <- read.table("weighted.field.unifrac", header = T, row.names = 1)
uuf.field.df <- read.table("unweighted.field.unifrac", header = T, row.names = 1)

### Weighted UniFrac
## Lundberg farms experiment
# Prinicipal coordinates analysis
wuf.field.pcoa <- pcoa(as.dist(wuf.field.df))
wuf.field.vectors <- cbind(field.map, wuf.field.pcoa$vectors)
head(wuf.field.pcoa$values)
wuf.field.vectors$site.label <- paste(wuf.field.vectors$Site, wuf.field.vectors$Cultivation, sep = " ")
wuf.field.vectors$comp.label <- paste(wuf.field.vectors$Compartment, wuf.field.vectors$Cultivation, sep = " ")

ggplot(wuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  theme_classic() +
  labs(x = "PCo1 (24.7%)", y = "PCo2 (16.9%)", legend = "", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) + 
  theme_classic() +
  labs(x = "PCo1 (24.7%)", y = "PCo2 (16.9%)", color = "Compartment") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = lat)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_gradientn(colours = rainbow(7)) + 
  theme_classic() +
  labs(x = "PCo1 (24.7%)", y = "PCo2 (16.9%)", color = "Latitude") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.field.vectors, aes(x = Axis.2, y = Axis.3, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_classic() +
  labs(x = "PCo2 (16.9%)", y = "PCo3 (10.7%)", legend = "", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.field.vectors, aes(x = Axis.1, y = Axis.3, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo1 (24.7%)", y = "PCo3 (10.7%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(wuf.field.vectors, aes(x = Axis.2, y = Axis.3, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo2 (16.9%)", y = "PCo3 (10.7%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.field.df) ~ Cultivation * Site * Compartment, data = field.map)

### Unweighted UniFrac
## Lundberg farms experiment
# Format map for experiment

# Prinicipal coordinates analysis
uuf.field.pcoa <- pcoa(as.dist(uuf.field.df))
uuf.field.vectors <- cbind(field.map, uuf.field.pcoa$vectors)
head(uuf.field.pcoa$values)
uuf.field.vectors$site.label <- paste(uuf.field.vectors$Site, uuf.field.vectors$Cultivation, sep = " ")
uuf.field.vectors$comp.label <- paste(uuf.field.vectors$Compartment, uuf.field.vectors$Cultivation, sep = " ")

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = site.label, shape = site.label)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_shape_manual(values = c(16, 16, 16, 16, 17, 16, 17, 17)) +
  theme_bw() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  #scale_shape_manual(values = c(16, 16, 16, 16, 17, 16, 17, 17)) +
  theme_classic() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

uuf.field.vectors$lat <- factor(uuf.field.vectors$lat)
ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = lat)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = brewer.pal(n = 10, "YlGnBu")[3:10]) + 
  theme_classic() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", color = "Latitude") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) + 
  #scale_shape_manual(values = c(17, 16)) +
  theme_classic() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "Compartment", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())


ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = lat, size = lat)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_gradientn(colours = rainbow(7)) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_classic() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "Latitude", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())
  

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo1 (11.8%)", y = "PCo2 (8.1%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.field.vectors, aes(x = Axis.1, y = Axis.3, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_color_manual(values = c("#377EB8", "#4DAF4A", "#984EA3")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo1 (11.8%)", y = "PCo3 (7.0%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.field.vectors, aes(x = Axis.2, y = Axis.3, color = Cultivation, shape = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  scale_shape_manual(values = c(17, 16)) +
  theme_bw() +
  labs(x = "PCo2 (8.1%)", y = "PCo3 (7.0%)", legend = "", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

uuf.field.vectors$Cultivation <- gsub("EcoFarm", "Ecofarm", uuf.field.vectors$Cultivation)
ggplot(uuf.field.vectors, aes(x = Axis.2, y = Axis.3, color = Cultivation)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = c("brown4", "darkorange")) + 
  theme_classic() +
  labs(x = "PCo2 (8.1%)", y = "PCo3 (7.0%)", legend = "", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())


x
## PERMUTATIONAL MANOVA
adonis(as.dist(uuf.field.df) ~ Cultivation * Site * Compartment, data = field.map)
