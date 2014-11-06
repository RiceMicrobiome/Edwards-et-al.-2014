library(vegan)
library(ggplot2)

## Load in the data
load("~/RMB/Publication/TimeCourse/TC+GH/Data/tc_gh_data.rda")
otu.table <- tc.gh.data[[1]]
map <- tc.gh.data[[2]]
tax <- tc.gh.data[[3]]

temp.map <- subset(map, Cultivation == "Greenhouse" & Site == "Davis" & Compartment == "Bulk Soil")
dav.map <- rbind(temp.map, subset(map, Site == "Davis" & Days != 21 & Days != 55 & Days != 34 & Cultivar == "M104"))
dav.map$Days[dav.map$Days == 30] = 42
dav.otu <- otu.table[,match(row.names(dav.map), colnames(otu.table))]
row.names(dav.otu) <- paste("Otu", row.names(dav.otu), sep = "")

## Get alpha diveristies
dav.adiv <- cbind(dav.map, Shannon = diversity(t(dav.otu)))
dav.adiv$Compartment <- gsub("Bulk_Soil", "Bulk Soil", dav.adiv$Compartment)
dav.adiv$Compartment <- factor(dav.adiv$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
dav.adiv <- subset(dav.adiv, Days <= 13)
summary(aov(exp(Shannon) ~ Compartment * Days, data = dav.adiv))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Bulk Soil")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizosphere")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizoplane")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Endosphere")))

kruskal.test(exp(Shannon) ~ Compartment * Days, data = dav.adiv)
kruskal.test(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Bulk Soil"))
kruskal.test(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizosphere"))
kruskal.test(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizoplane"))
s
## Plot the results
comp.cols <- c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")
ggplot(subset(dav.adiv, Days > 0), aes(x = "", y = exp(Shannon), fill = Compartment)) +
  geom_boxplot(width = 1) +
  scale_fill_manual(values = comp.cols) +
  facet_grid(Compartment ~ Days) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(text = element_text(size = 20), axis.ticks.x = element_blank(), legend.key = element_blank())

ggplot(subset(dav.adiv, Days > 0), aes(x = factor(Days), y = exp(Shannon), fill = Compartment)) +
  geom_boxplot(width = 1) +
  scale_fill_manual(values = comp.cols) +
  facet_grid(Compartment ~ .) +
  theme_bw() +
  labs(x = "", y = "Effective Species") +
  theme(text = element_text(size = 20), legend.key = element_blank())

summary(aov(exp(Shannon) ~ Compartment * Days, data = subset(dav.adiv, Days > 0 & Days < 42)))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Endosphere" & Days > 0 & Days < 42)))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizoplane" & Days > 0 & Days < 42)))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizosphere" & Days > 0 & Days < 42)))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Bulk Soil" & Days > 0 & Days < 42)))

summary(aov(exp(Shannon) ~ Compartment * Days, data = subset(dav.adiv, Days > 0 & Cultivation == "TC")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Endosphere" & Days > 0 & Cultivation =="TC")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizoplane" & Days > 0 & Cultivation =="TC")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Rhizosphere" & Days > 0 & Cultivation =="TC")))
summary(aov(exp(Shannon) ~ Days, data = subset(dav.adiv, Compartment == "Bulk Soil" & Days > 0 & Cultivation =="TC")))
