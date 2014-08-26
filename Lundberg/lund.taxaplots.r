library(reshape)
library(ggplot2)

### Set the working directory
setwd("~/RICE/Publication/Data/FieldExp/")

### Load in the data
field.counts <- read.table("field_otu_table.txt", header = T, row.names = 1)
field.map <- read.table("field_map.txt", header = T, row.names = 1)
field.tax <- read.table("field_tax.txt", header = T, row.names = 1)

### Format the map a bit
field.map$Field <- NULL
field.map$Run <- NULL

############################################
## Phylum Plots
############################################
### Aggregate taxa and counts based on phyla
field.phy <- aggregate(field.counts, data.frame(field.tax$Phylum), sum)
row.names(field.phy) <- field.phy[,1]
field.phy <- field.phy[,-1]

### Extract the top15 highest represented phyla
top.15 <- names(head(sort(rowSums(field.phy), decreasing = T), 15))
other <- colSums(field.phy[!row.names(field.phy)%in%top.15,])
field.phy.15 <- rbind(field.phy[row.names(field.phy)%in%top.15,], other = other)

### Melt into a long data frame for plotting
field.phy.whole <- melt(cbind(field.map, t(field.phy.15)))
field.phy.whole$Compartment <- factor(field.phy.whole$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))

### Define the colors
colors <- c("mediumpurple3", "palegreen4", "skyblue1","darkorange", "blue","darkseagreen", "chocolate4","rosybrown1", "red", "gold", "dodgerblue4", "orange", "forestgreen", "grey", "orchid",  "grey50")

### Plot this thing
ggplot(field.phy.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 30), axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 40, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))

### Plot Proteobacterial classes
prot.tax <- subset(field.tax, Phylum == "Proteobacteria")
prot.counts <- field.counts[match(row.names(prot.tax), row.names(field.counts)),]
field.prot <- aggregate(prot.counts, data.frame(prot.tax$Class), sum)
row.names(field.prot) <- field.prot[,1]
field.prot <- field.prot[,-1]
field.prot.whole <- melt(cbind(field.map, t(field.prot)))

prot.cols <- c(colors[1:4], colors[6], "grey")
field.prot.whole$Compartment <- factor(field.prot.whole$Compartment, levels = c("Rhizosphere", "Rhizoplane", "Endosphere"))
ggplot(field.prot.whole, aes(x = Compartment, y = value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_grid(.~ Site) +
  scale_fill_manual(values = prot.cols) +
  labs(x = "", y = "Proportion of Counts", fill = "Phylum") +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 10), axis.text.x = element_text(angle = 30, hjust = 1)) +
  guides(fill = guide_legend(reverse = T))