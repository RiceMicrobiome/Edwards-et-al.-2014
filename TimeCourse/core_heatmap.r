library(ggplot2)
library(reshape2)

source("~/RMB/Publication/Scripts/General/rmb_functions.r")
setwd("~/RMB/Publication/TimeCourse/TC+GH/Data/")

## Load in relevant files
load("tc_gh_data.rda")
counts <- tc.gh.data[[1]]
map <- tc.gh.data[[2]]
core.otus <- read.table("endo_core_gh.txt", header = T)

## Format everything
good.map <- subset(map, Days <= 13)
good.map$Compartment <- gsub("Bulk_Soil", "Bulk Soil", good.map$Compartment)
row.names(counts) <- paste("Otu", row.names(counts), sep = "")
good.counts <- counts[match(core.otus$e.int, row.names(counts)), match(row.names(good.map), colnames(counts))]
good.counts <- good.counts[apply(good.counts, 1, sum) > 30,]

clusts <- hclust(dist(1-cor(t(good.counts), method = "spearman")), method = "average")
clusts.order <- data.frame(cbind(OTU = clusts$labels, Order = clusts$order))
clusts.order$Order <- as.numeric(as.character(clusts.order$Order))
target <- 1:max(clusts.order$Order)
clusts.order <- clusts.order[match(target, clusts.order$Order),]

## Get mean of sample types
good.map$sampType <- with(good.map, paste(Compartment, Days, sep = "_"))
type.counts <- aggregate(t(good.counts), data.frame(good.map$sampType), mean)
row.names(type.counts) <- type.counts[,1]
type.counts <- type.counts[,-1]
type.counts <- scale(type.counts)
whole.data <- cbind(colsplit(row.names(type.counts), pattern = "_", names = c("Compartment", "Days")), type.counts)
whole.data$Days <- factor(whole.data$Days, levels = c(0,1,2,3,5,8,13))
whole.data$Compartment <- factor(whole.data$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
whole.data.m <- melt(whole.data)
whole.data.m$variable <- factor(whole.data.m$variable, levels = as.character(clusts.order$OTU))

## Plot it
ggplot(whole.data.m, aes(x = Days, y = variable, fill = value)) +
  geom_tile() +
  facet_grid(Compartment ~.) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  scale_fill_gradient2(low = 'black', mid = 'yellow', high = 'orange',midpoint = 2)

## Endosphere alone
endo.data <- subset(whole.data, Compartment == "Endosphere")
endo.data$Days <- factor(endo.data$Days, levels = c(1,2,3,5,8,13))
endo.data.m <- melt(endo.data)
ggplot(endo.data.m, aes(x = Days, y = variable, fill = log2(value + 1))) +
  geom_tile()

## Rhizoplane alone
rp.data <- subset(whole.data, Compartment == "Rhizoplane")
rp.data$Days <- factor(rp.data$Days, levels = c(1,2,3,5,8,13))
rp.data.m <- melt(rp.data)
ggplot(rp.data.m, aes(x = Days, y = variable, fill = log2(value + 1))) +
  geom_tile()

## Rhizosphere alone
rs.data <- subset(whole.data, Compartment == "Rhizosphere")
rs.data$Days <- factor(rs.data$Days, levels = c(1,2,3,5,8,13))
rs.data.m <- melt(rs.data)
ggplot(rs.data.m, aes(x = Days, y = variable, fill = log2(value + 1))) +
  geom_tile()
