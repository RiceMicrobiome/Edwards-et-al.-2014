library(ggplot2)
library(vegan)
library(ape)

setwd("~/RMB/TimeCourse/TC+GH/Data/")

wuf <- read.table("weighted_unifrac_GH_TC_norm.txt", header = T, row.names = 1)
uuf <- read.table("unweighted_unifrac_GH_TC_norm.txt", header = T, row.names = 1)
map <- read.table("TC+GH.map", header = T, row.names = 1, sep = "\t")
map <- map[match(row.names(wuf), row.names(map)),]
map$Compartment <- gsub("Bulk_Soil", "Bulk Soil", map$Compartment)
map$Days[which(map$Cultivation == "Greenhouse")] <- 42

wuf.pc <- pcoa(as.dist(wuf))
head(wuf.pc$values)
wuf.data <- cbind(wuf.pc$vectors[,1:4], map)
wuf.data$Compartment <- factor(wuf.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
comp.cols <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")

ggplot(wuf.data, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 9, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (32.9%)", y = "PCo 2 (15.8%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(wuf.data, aes(x = Axis.1, y = Axis.2, color = Days, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_classic() +
  labs(x = "PCo 1 (32.9%)", y = "PCo 2 (15.8%)") +
  scale_color_gradientn(colours = rainbow(10)) +
  theme(text = element_text(size = 20), legend.key = element_blank())

uuf.pc <- pcoa(as.dist(uuf))
head(uuf.pc$values)
uuf.data <- cbind(uuf.pc$vectors, map)
uuf.data$Compartment <- factor(uuf.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))

ggplot(uuf.data, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (15.4%)", y = "PCo 2 (12.8%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.data, aes(x = Axis.1, y = Axis.2, color = Days, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (15.4%)", y = "PCo 2 (12.8%)") +
  scale_color_gradientn(colours = rainbow(10)) +
  theme(text = element_text(size = 20), legend.key = element_blank())


pcoa_sub <- function(distance, map) {
  distance <- distance[match(row.names(map), row.names(distance)), match(row.names(map), colnames(distance))]
  pc <- pcoa(distance)
  return(pc)
}

dav.wpc <- pcoa_sub(wuf, subset(map, Site == "Davis"))
dav.wdata <- cbind(subset(map, Site == "Davis"), dav.wpc$vectors[,1:4])
head(dav.wpc$values)
dav.wdata$Compartment <- factor(dav.wdata$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(dav.wdata, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (36.0%)", y = "PCo 2 (14.4%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(dav.wdata, aes(x = Axis.1, y = Axis.2, color = Days)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (36.0%)", y = "PCo 2 (14.4%)") +
  scale_color_gradientn(colours = rainbow(10)) +
  theme(text = element_text(size = 20), legend.key = element_blank())

dav.upc <- pcoa_sub(uuf, subset(map, Site == "Davis"))
dav.udata <- cbind(subset(map, Site == "Davis"), dav.upc$vectors[,1:4])
head(dav.upc$values)
dav.udata$Compartment <- factor(dav.udata$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(dav.udata, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (21.4%)", y = "PCo 2 (10.0%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(dav.udata, aes(x = Axis.1, y = Axis.2, color = Days)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (21.4%)", y = "PCo 2 (10.0%)") +
  scale_color_gradientn(colours = rainbow(10)) +
  theme(text = element_text(size = 20), legend.key = element_blank())

### Remove timepoint 21 and above
map.set <- subset(map, Days != 21 & Days != 34 & Days != 55)
wuf.set <- wuf[match(row.names(map.set), row.names(wuf)), match(row.names(map.set), colnames(wuf))]
uuf.set <- uuf[match(row.names(map.set), row.names(uuf)), match(row.names(map.set), colnames(uuf))]

## Weighted
wuf.set.pc <- pcoa(as.dist(wuf.set))
wuf.set.data <- cbind(map.set, wuf.set.pc$vectors[,1:4])
wuf.set.data$Compartment <- factor(wuf.set.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
head(wuf.set.pc$values)

ggplot(wuf.set.data, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (34.1%)", y = "PCo 2 (16.2%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(wuf.set.data, aes(x = Axis.1, y = Axis.2, color = Days, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (34.1%)", y = "PCo 2 (16.2%)") +
  scale_color_gradientn(colours = rainbow(7)) +
  theme(text = element_text(size = 20), legend.key = element_blank())

uuf.set.pc <- pcoa(as.dist(uuf.set))
uuf.set.data <- cbind(map.set, uuf.set.pc$vectors[,1:4])
uuf.set.data$Compartment <- factor(uuf.set.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
head(uuf.set.pc$values)
ggplot(uuf.set.data, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (15.9%)", y = "PCo 2 (12.9%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

### Only do days greater than 34
lower.map <- subset(map, Days != 21 &Days != 34 & Days != 55)
lower.wuf <- wuf[match(row.names(lower.map), row.names(wuf)), match(row.names(lower.map), colnames(wuf))]
lower.uuf <- uuf[match(row.names(lower.map), row.names(uuf)), match(row.names(lower.map), colnames(uuf))]

lower.wuf.pc <- pcoa(as.dist(lower.wuf))
lower.wuf.data <- cbind(lower.map, lower.wuf.pc$vectors[,1:4])
lower.wuf.data$Compartment <- factor(lower.wuf.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
head(lower.wuf.pc$values)
ggplot(lower.wuf.data, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "PCo 1 (34.1%)", y = "PCo 2 (16.2%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(lower.wuf.data, aes(x = Axis.1, y = Axis.2, color = Days, shape = Site)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 9, alpha = 0.9) +
  theme_classic() +
  labs(x = "PCo 1 (37.3%)", y = "PCo 2 (16.9%)") +
  scale_color_gradientn(colours = rainbow(3)) +
  theme(text = element_text(size = 30), legend.key = element_blank())


### Only do Davis greater than 34
dav.map <- subset(map, Site == "Davis" & Days != 21 &Days != 34 & Days != 55)
dav.wuf <- wuf[match(row.names(dav.map), row.names(wuf)), match(row.names(dav.map), colnames(wuf))]
dav.uuf <- uuf[match(row.names(dav.map), row.names(uuf)), match(row.names(dav.map), colnames(uuf))]

dav.wuf.pc <- pcoa(as.dist(dav.wuf))
dav.wuf.data <- cbind(dav.map, dav.wuf.pc$vectors[,1:4])
dav.wuf.data$Compartment <- factor(dav.wuf.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
head(dav.wuf.pc$values)
ggplot(dav.wuf.data, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 9, alpha = 0.9, shape = 17) +
  theme_classic() +
  labs(x = "PCo 1 (37.3%)", y = "PCo 2 (16.9%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(dav.wuf.data, aes(x = Axis.1, y = Axis.2, color = Days)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 9, alpha = 0.9, shape = 17) +
  theme_classic() +
  labs(x = "PCo 1 (37.3%)", y = "PCo 2 (16.9%)") +
  scale_color_gradientn(colours = rainbow(5), guide = guide_legend(direction = "horizontal", label.position = "bottom", title.position = "top")) +
  theme(text = element_text(size = 30), legend.key = element_blank())

dav.uuf.pc <- pcoa(as.dist(dav.uuf))
dav.uuf.data <- cbind(dav.map, dav.uuf.pc$vectors[,1:4])
dav.uuf.data$Compartment <- factor(dav.uuf.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
head(dav.uuf.pc$values)
ggplot(dav.uuf.data, aes(x = Axis.1, y = Axis.2, color = Compartment)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (22.4%)", y = "PCo 2 (12.2%)") +
  scale_color_manual(values = comp.cols) +
  theme(text = element_text(size = 20), legend.key = element_blank())

ggplot(dav.uuf.data, aes(x = Axis.1, y = Axis.2, color = Days)) +
  geom_vline(x = 0, alpha = 0.5) +
  geom_hline(y = 0, alpha = 0.5) +
  geom_point(size = 7, alpha = 0.9) +
  theme_bw() +
  labs(x = "PCo 1 (22.4%)", y = "PCo 2 (12.2%)") +
  scale_color_gradientn(colours = rainbow(5)) +
  theme(text = element_text(size = 20), legend.key = element_blank())
