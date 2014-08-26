library(ggplot2)
library(vegan)

setwd("~/RMB/Publication/Data/FieldExp/")

counts <- read.table("field_otu_table.txt", header = T, row.names = 1)
map <- read.table("field_map.txt", header = T, row.names = 1)
map <- map[match(colnames(counts), row.names(map)),]

adiv <- data.frame(cbind(map, Shannon = diversity(t(counts))))
summary(aov(exp(Shannon) ~ Cultivation * Compartment * Site, data = adiv))

t.test(subset(adiv, Compartment == "Rhizosphere" & Cultivation == "Organic")$Shannon,
       subset(adiv, Compartment == "Rhizosphere" & Cultivation == "EcoFarm")$Shannon)

t.test(subset(adiv, Compartment == "Rhizoplane" & Cultivation == "Organic")$Shannon,
       subset(adiv, Compartment == "Rhizoplane" & Cultivation == "EcoFarm")$Shannon)

t.test(subset(adiv, Compartment == "Endosphere" & Cultivation == "Organic")$Shannon,
       subset(adiv, Compartment == "Endosphere" & Cultivation == "EcoFarm")$Shannon)

counts <- read.table("../GreenhouseExp/gh_otu_table.txt", header = T, row.names = 1)
map <- read.table ("../GreenhouseExp/gh_map.txt", header = T, row.names = 1)
summary(aov(exp(Shannon) ~ Compartment * Cultivar * Site, data = adiv))
summary(aov(exp(Shannon) ~ Compartment * Cultivar * Site, data = subset(adiv, Compartment != "Bulk Soil")))
