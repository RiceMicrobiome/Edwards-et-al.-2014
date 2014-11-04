library(ggplot2)

## Load the data
load("~/RMB/Publication/Data/GreenhouseExp/comp.otus.rda")
comp.otus$logFC <- comp.otus$logFC * -1
tax <- read.table("~/RMB/Publication/Data/GreenhouseExp/gh_tax.txt", header = T, row.names = 1)

## Find which OTUs are significant
comp.otus$sig <- ifelse(comp.otus$padj <= 0.01, "sig", "ns")

## Put phyla to plot in
int.phyla <- c("Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Chloroflexi", "Fibrobacteres", "Firmicutes", "Gemmatimonadetes", "Nitrospirae", "Planctomycetes", "Proteobacteria", "Spirochaetes", "unclassified", "Verrucomicrobia", "WS3")
comp.otus$Phylum <- tax[match(comp.otus$OTU, tax$OTU),]$Phylum
levels(comp.otus$Phylum) <- c(levels(comp.otus$Phylum), "Other", "ns")
comp.otus$Phylum[!comp.otus$Phylum%in%int.phyla] <- "Other"
comp.otus$Phylum[which(comp.otus$sig == "ns")] <- "ns"

## Get everything straight for plotting
phy.cols <- c("mediumpurple3", "palegreen4", "skyblue1", "darkorange", "darkseagreen", "firebrick4", "red", "gold", "goldenrod4", "dodgerblue4", "orange", "forestgreen", "grey", "orchid", "darkmagenta", "black", "grey50")

## Plot it
ggplot(comp.otus, aes(x = logCPM, y = logFC, color = Phylum, alpha = sig)) +
  geom_point() +
  scale_color_manual(values = phy.cols) +
  facet_grid(.~Comp) +
  theme_bw() +
  theme(text = element_text(size = 20))