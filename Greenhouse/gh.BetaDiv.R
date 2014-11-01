library(vegan)
library(reshape2)
library(ggplot2)
library(ape)


### Set working directory for this analysis
setwd("~/RMB/Publication/Data/GreenhouseUniFrac/")

### Load in the data
# Load in the mapping file and format it a bit
gh.map <- read.table("gh_map.txt",header = T, row.names = 1)

# Load in weighted and unweighted unifrac distances
wuf.gh.df <- read.table("weighted.gh.unifrac", header = T, row.names = 1)
uuf.gh.df <- read.table("unweighted.gh.unifrac", header = T, row.names = 1)

##########################################################
## Weighted UniFrac
##########################################################

# Pull out relevant data in the weighted distance matrix
wuf.gh.df <- wuf.gh.df[match(row.names(gh.map), row.names(wuf.gh.df)), match(row.names(gh.map), colnames(wuf.gh.df))]

# Prinicipal coordinates analysis
wuf.gh.pcoa <- pcoa(as.dist(wuf.gh.df))
wuf.gh.vectors <- cbind(gh.map, wuf.gh.pcoa$vectors)
head(wuf.gh.pcoa$values)
wuf.gh.vectors$site.comp <- paste(wuf.gh.vectors$Compartment, wuf.gh.vectors$Site, sep = " ")

wuf.gh.vectors$Compartment <- factor(wuf.gh.vectors$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
wuf.gh.vectors$Cultivar <- gsub("Nipponebare", "Nipponbare", wuf.gh.vectors$Cultivar)
wuf.gh.vectors$Cultivar <- factor(wuf.gh.vectors$Cultivar, levels = c("Soil", "Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))

comp.cols <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")
cult.cols <- c("black", "#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")

ggplot(wuf.gh.vectors, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = comp.cols) +
  theme_classic() +
  labs(x = "PCo1 (46.3%)", y = "PCo2 (11.5%)", color = "Rhizocompartment", shape = "Site") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.gh.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivar, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = cult.cols) +
  theme_classic() +
  labs(x = "PCo1 (46.3%)", y = "PCo2 (11.5%)", color = "Rhizocompartment", shape = "Site") +
  theme(text = element_text(size = 30), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.gh.df) ~ Site * Cultivar * Compartment + Field, strata = gh.map$Field, data = gh.map)
adonis(as.dist(wuf.gh.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(gh.map$Site, gh.map$Field), data = gh.map)

### Only the rhizosphere compartment vs bulk soil
wuf.gh.rs.map <- subset(gh.map, Compartment == "Rhizosphere" | Compartment == "Bulk Soil")
wuf.gh.rs.df <- wuf.gh.df[match(row.names(wuf.gh.rs.map), row.names(wuf.gh.df)), match(row.names(wuf.gh.rs.map), colnames(wuf.gh.df))]

wuf.gh.rs.map <- subset(gh.map, Compartment == "Rhizosphere")
wuf.gh.rs.df <- wuf.gh.df[match(row.names(wuf.gh.rs.map), row.names(wuf.gh.df)), match(row.names(wuf.gh.rs.map), colnames(wuf.gh.df))]

# Principal coordinate analysis
wuf.gh.rs.pcoa <- pcoa(as.dist(wuf.gh.rs.df))
wuf.gh.rs.vectors <- cbind(wuf.gh.rs.map, wuf.gh.rs.pcoa$vectors)
head(wuf.gh.rs.pcoa$values)
wuf.gh.rs.vectors$site.comp <- paste(wuf.gh.rs.vectors$Site, wuf.gh.rs.vectors$Compartment, sep = " ")
wuf.gh.rs.vectors$site.cult <- factor(paste(wuf.gh.rs.vectors$Site, wuf.gh.rs.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                       "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                       "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))                                    

cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")
cult.cols <- c("#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")
cult.cols <- c("#006600", "#B2D1B2", "#FF9933", "#FFD6AD", "#660066", "#D1B2D1")

ggplot(wuf.gh.rs.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#984EA3"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (27.0%)", y = "PCo2 (20.8%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

wuf.gh.rs.vectors$Cultivar <- gsub("Glab_B", "Glab B", wuf.gh.rs.vectors$Cultivar)
wuf.gh.rs.vectors$Cultivar <- gsub("Glab_E", "Glab E", wuf.gh.rs.vectors$Cultivar)
wuf.gh.rs.vectors$Cultivar <- gsub("Nipponebare", "Nipponbare", wuf.gh.rs.vectors$Cultivar)
wuf.gh.rs.vectors$Cultivar <- factor(wuf.gh.rs.vectors$Cultivar, levels = c("Glab B", "Glab E", "93-11", "IR50", "M104", "Nipponbare"))

ggplot(wuf.gh.rs.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivar, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  #scale_shape_manual(values = c(rep(16, 6), rep(17, 6), rep(15, 6))) +
  theme_classic() +
  labs(x = "PCo1 (30.4%)", y = "PCo2 (24.2%)", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.gh.rs.vectors, aes(x = Axis.2, y = Axis.3, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 6), rep(17, 6), rep(15, 6))) +
  theme_bw() +
  labs(x = "PCo2 (24.3%)", y = "PCo3 (18.5%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.gh.rs.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = wuf.gh.rs.map$Site, data = wuf.gh.rs.map)

### Only the rhizoplane compartment
wuf.gh.rp.map <- subset(gh.map, Compartment == "Rhizoplane" | Compartment == "Bulk Soil")
wuf.gh.rp.df <- wuf.gh.df[match(row.names(wuf.gh.rp.map), row.names(wuf.gh.df)), match(row.names(wuf.gh.rp.map), colnames(wuf.gh.df))]

wuf.gh.rp.map <- subset(gh.map, Compartment == "Rhizoplane")
wuf.gh.rp.df <- wuf.gh.df[match(row.names(wuf.gh.rp.map), row.names(wuf.gh.df)), match(row.names(wuf.gh.rp.map), colnames(wuf.gh.df))]


# Principal coordinate analysis
wuf.gh.rp.pcoa <- pcoa(as.dist(wuf.gh.rp.df))
wuf.gh.rp.vectors <- cbind(wuf.gh.rp.map, wuf.gh.rp.pcoa$vectors)
head(wuf.gh.rp.pcoa$values)
wuf.gh.rp.vectors$site.comp <- paste(wuf.gh.rp.vectors$Site, wuf.gh.rp.vectors$Compartment, sep = " ")
wuf.gh.rp.vectors$site.cult <- factor(paste(wuf.gh.rp.vectors$Site, wuf.gh.rp.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                   "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                   "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))


cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")
cult.cols <- c("#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")

ggplot(wuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#4DAF4A"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (35.9%)", y = "PCo2 (15.5%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())
wuf.gh.rp.vectors$Cultivar <- gsub("Nipponebare", "Nipponbare", wuf.gh.rp.vectors$Cultivar)
wuf.gh.rp.vectors$Cultivar <- factor(wuf.gh.rp.vectors$Cultivar, levels = c("Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))
ggplot(wuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivar, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = cult.cols) +
  theme_classic() +
  labs(x = "PCo1 (35.9%)", y = "PCo2 (15.5%)", color = "Cultivar", shape = "Soil Source") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.3, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (28.1%)", y = "PCo3 (25.2%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.gh.rp.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(wuf.gh.rp.map$Site, wuf.gh.rp.map$Field), data = wuf.gh.rp.map)


### Only the endosphere compartment
wuf.gh.e.map <- subset(gh.map, Compartment == "Endosphere")
wuf.gh.e.df <- wuf.gh.df[match(row.names(wuf.gh.e.map), row.names(wuf.gh.df)), match(row.names(wuf.gh.e.map), colnames(wuf.gh.df))]

# Principal coordinate analysis
wuf.gh.e.pcoa <- pcoa(as.dist(wuf.gh.e.df))
wuf.gh.e.vectors <- cbind(wuf.gh.e.map, wuf.gh.e.pcoa$vectors)
head(wuf.gh.e.pcoa$values)
wuf.gh.e.vectors$site.comp <- paste(wuf.gh.e.vectors$Site, wuf.gh.e.vectors$Compartment, sep = " ")
wuf.gh.e.vectors$site.cult <- factor(paste(wuf.gh.e.vectors$Site, wuf.gh.e.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                   "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                   "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))


cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")

ggplot(wuf.gh.e.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#377EB8"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (55.2%)", y = "PCo2 (10.3%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

wuf.gh.e.vectors$Cultivar <- gsub("Nipponebare", "Nipponbare", wuf.gh.e.vectors$Cultivar)
wuf.gh.e.vectors$Cultivar <- factor(wuf.gh.e.vectors$Cultivar, levels = c("Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))
ggplot(wuf.gh.e.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivar, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = cult.cols) +
  theme_classic() +
  labs(x = "PCo1 (29.0%)", y = "PCo2 (21.9%)", color = "Cultivar", shape = "Soil Source") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(wuf.gh.e.vectors, aes(x = Axis.2, y = Axis.3, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (55.2%)", y = "PCo3 (7.0%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.gh.e.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(wuf.gh.e.map$Site, wuf.gh.e.map$Field), data = wuf.gh.e.map)


##########################################################
## Unweighted UniFrac
##########################################################

# Pull out relevant data in the weighted distance matrix
uuf.gh.df <- uuf.gh.df[match(row.names(gh.map), row.names(uuf.gh.df)), match(row.names(gh.map), colnames(uuf.gh.df))]

# Prinicipal coordinates analysis
uuf.gh.pcoa <- pcoa(as.dist(uuf.gh.df))
uuf.gh.vectors <- cbind(gh.map, uuf.gh.pcoa$vectors)
head(uuf.gh.pcoa$values)
uuf.gh.vectors$site.comp <- paste(uuf.gh.vectors$Site, uuf.gh.vectors$Compartment, sep = " ")
uuf.gh.vectors$Compartment <- factor(uuf.gh.vectors$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
uuf.gh.vectors$Cultivar <- gsub("Nipponebare", "Nipponbare", uuf.gh.vectors$Cultivar)
uuf.gh.vectors$Cultivar <- factor(uuf.gh.vectors$Cultivar, levels = c("Soil", "Glab_B", "Glab_E", "93-11", "IR50", "M104", "Nipponbare"))

ggplot(uuf.gh.vectors, aes(x = Axis.1, y = Axis.2, color = Compartment, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A",  "#984EA3")) +
  theme_classic() +
  labs(x = "PCo1 (18.1%)", y = "PCo2 (14.9%)", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

ggplot(uuf.gh.vectors, aes(x = Axis.1, y = Axis.2, color = Cultivar, shape = Site)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 9, alpha= 0.75) +
  scale_color_manual(values = cult.cols) +
  theme_classic() +
  labs(x = "PCo1 (18.1%)", y = "PCo2 (14.9%)", color = "", shape = "") +
  theme(text = element_text(size = 30), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(wuf.gh.df) ~ Site * Cultivar * Compartment + Field, strata = gh.map$Field, data = gh.map)
adonis(as.dist(uuf.gh.df) ~ Compartment * Site * Cultivar + paste(gh.map$Site, gh.map$Field), strata = paste(gh.map$Site, gh.map$Field), data = gh.map)

### Only the rhizosphere compartment vs bulk soil
uuf.gh.rs.map <- subset(gh.map, Compartment == "Rhizosphere" | Compartment == "Bulk Soil")
uuf.gh.rs.df <- uuf.df[match(row.names(uuf.gh.rs.map), row.names(uuf.gh.df)), match(row.names(uuf.gh.rs.map), colnames(uuf.gh.df))]

# Principal coordinate analysis
uuf.gh.rs.pcoa <- pcoa(as.dist(uuf.gh.rs.df))
uuf.gh.rs.vectors <- cbind(uuf.gh.rs.map, uuf.gh.rs.pcoa$vectors)
head(uuf.gh.rs.pcoa$values)
uuf.gh.rs.vectors$site.comp <- paste(uuf.gh.rs.vectors$Site, uuf.gh.rs.vectors$Compartment, sep = " ")
uuf.gh.rs.vectors$site.cult <- factor(paste(uuf.gh.rs.vectors$Site, uuf.gh.rs.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                       "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                       "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))                                    

cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")

ggplot(uuf.gh.rs.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#984EA3"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (27.0%)", y = "PCo2 (12.1%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.rs.vectors, aes(x = Axis.1, y = Axis.2, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (27.0%)", y = "PCo2 (12.1%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.rs.vectors, aes(x = Axis.1, y = Axis.3, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (27.0%)", y = "PCo3 (5.5%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(uuf.gh.rs.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(uuf.gh.rs.map$Site, uuf.gh.rs.map$Field), data = uuf.gh.rs.map)

### Only the rhizoplane compartment
uuf.gh.rp.map <- subset(gh.map, Compartment == "Rhizoplane" | Compartment == "Bulk Soil")
uuf.gh.rp.df <- uuf.df[match(row.names(uuf.gh.rp.map), row.names(uuf.gh.df)), match(row.names(uuf.gh.rp.map), colnames(uuf.gh.df))]

# Principal coordinate analysis
uuf.gh.rp.pcoa <- pcoa(as.dist(uuf.gh.rp.df))
uuf.gh.rp.vectors <- cbind(uuf.gh.rp.map, uuf.gh.rp.pcoa$vectors)
head(uuf.gh.rp.pcoa$values)
uuf.gh.rp.vectors$site.comp <- paste(uuf.gh.rp.vectors$Site, uuf.gh.rp.vectors$Compartment, sep = " ")
uuf.gh.rp.vectors$site.cult <- factor(paste(uuf.gh.rp.vectors$Site, uuf.gh.rp.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                       "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                       "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))


cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")

ggplot(uuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#4DAF4A"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (19.1%)", y = "PCo2 (9.54%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.2, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (19.1%)", y = "PCo2 (9.54%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.rp.vectors, aes(x = Axis.1, y = Axis.4, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (19.1%)", y = "PCo3 (7.2%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(uuf.gh.rp.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(uuf.gh.rp.map$Site, uuf.gh.rp.map$Field), data = uuf.gh.rp.map)

### Only the endosphere compartment
uuf.gh.e.map <- subset(gh.map, Compartment == "Endosphere" | Compartment == "Bulk Soil")
uuf.gh.e.df <- uuf.df[match(row.names(uuf.gh.e.map), row.names(uuf.gh.df)), match(row.names(uuf.gh.e.map), colnames(uuf.gh.df))]

# Principal coordinate analysis
uuf.gh.e.pcoa <- pcoa(as.dist(uuf.gh.e.df))
uuf.gh.e.vectors <- cbind(uuf.gh.e.map, uuf.gh.e.pcoa$vectors)
head(uuf.gh.e.pcoa$values)
uuf.gh.e.vectors$site.comp <- paste(uuf.gh.e.vectors$Site, uuf.gh.e.vectors$Compartment, sep = " ")
uuf.gh.e.vectors$site.cult <- factor(paste(uuf.gh.e.vectors$Site, uuf.gh.e.vectors$Cultivar, sep = " "), levels = c("Arbuckle Soil", "Arbuckle Glab_B",  "Arbuckle Glab_E", "Arbuckle IR50", "Arbuckle 93-11", "Arbuckle Nipponebare", "Arbuckle M104",
                                                                                                                    "Davis Soil", "Davis Glab_B", "Davis Glab_E", "Davis IR50", "Davis 93-11", "Davis Nipponebare", "Davis M104",
                                                                                                                    "Sacramento Soil", "Sacramento Glab_B", "Sacramento Glab_E", "Sacramento IR50", "Sacramento 93-11", "Sacramento Nipponebare", "Sacramento M104"))


cult.cols <- c("grey50", "#F56600", "#E8A779", "#008FF5", "#79BAE8", "#BA00BA", "#E87DE8")

ggplot(uuf.gh.e.vectors, aes(x = Axis.1, y = Axis.2, color = site.comp, shape = site.comp)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(c("#E41A1C", "#377EB8"), 3)) +
  scale_shape_manual(values = c(rep(16, 2), rep(17, 2), rep(15, 2))) +
  theme_bw() +
  labs(x = "PCo1 (21.7%)", y = "PCo2 (12.1%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.e.vectors, aes(x = Axis.1, y = Axis.2, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (21.7%)", y = "PCo2 (12.1%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

ggplot(uuf.gh.e.vectors, aes(x = Axis.1, y = Axis.3, color = site.cult, shape = site.cult)) +
  geom_hline(y = 0, alpha = 0.3) +
  geom_vline(x = 0, alpha = 0.3) +
  geom_point(size = 7.5, alpha= 0.75) +
  scale_color_manual(values = rep(cult.cols, 3)) +
  scale_shape_manual(values = c(rep(16, 7), rep(17, 7), rep(15, 7))) +
  theme_bw() +
  labs(x = "PCo1 (21.7%)", y = "PCo3 (6.1%)", color = "", shape = "") +
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.key = element_blank())

## PERMUTATIONAL MANOVA
adonis(as.dist(uuf.gh.e.df) ~ Compartment * Site * Cultivar + paste(Site, Field), strata = paste(uuf.gh.e.map$Site, uuf.gh.e.map$Field), data = uuf.gh.e.map)

#######################################
## Do analyses on different soils
#######################################
### Weighted
## Arbuckle
arb.map <- subset(gh.map, Site == "Arbuckle")
arb.wuf.dist <- wuf.gh.df[match(row.names(arb.map), row.names(wuf.gh.df)), match(row.names(arb.map), colnames(wuf.gh.df))]
adonis(as.dist(arb.wuf.dist) ~ Compartment * Cultivar + Field, data = arb.map)

## Sacramento
sac.map <- subset(gh.map, Site == "Sacramento")
sac.wuf.dist <- wuf.gh.df[match(row.names(sac.map), row.names(wuf.gh.df)), match(row.names(sac.map), colnames(wuf.gh.df))]
adonis(as.dist(sac.wuf.dist) ~ Compartment * Cultivar + Field, data = sac.map)

## Davis
dav.map <- subset(gh.map, Site == "Davis")
dav.wuf.dist <- wuf.gh.df[match(row.names(dav.map), row.names(wuf.gh.df)), match(row.names(dav.map), colnames(wuf.gh.df))]
adonis(as.dist(dav.wuf.dist) ~ Compartment * Cultivar + Field, data = dav.map)

### Unweighted
## Arbuckle
arb.uuf.dist <- uuf.gh.df[match(row.names(arb.map), row.names(uuf.gh.df)), match(row.names(arb.map), colnames(uuf.gh.df))]
adonis(as.dist(arb.uuf.dist) ~ Compartment * Cultivar + Field, data = arb.map)

## Sacramento
sac.uuf.dist <- uuf.gh.df[match(row.names(sac.map), row.names(uuf.gh.df)), match(row.names(sac.map), colnames(uuf.gh.df))]
adonis(as.dist(sac.uuf.dist) ~ Compartment * Cultivar + Field, data = sac.map)

## Davis
dav.uuf.dist <- uuf.gh.df[match(row.names(dav.map), row.names(uuf.gh.df)), match(row.names(dav.map), colnames(uuf.gh.df))]
adonis(as.dist(dav.uuf.dist) ~ Compartment * Cultivar + Field, data = dav.map)


