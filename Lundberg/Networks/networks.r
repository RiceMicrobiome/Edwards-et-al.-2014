library(dynamicTreeCut)

load("~/RMB/Publication/Data/FieldExp/interesting_counts.rda")
tax <- read.table("RMB/Publication/Data/FieldExp/field_tax.txt", header = T, row.names = 1)

## Get spearman correlations for interesting OTUs
otu.cor <- cor(lund.interesting, method = "pearson")

## Agglomerative clustering on OTUs
otu.tre <- hclust(as.dist(1 - otu.cor), method = "average")

## Dynamic pruning of tree
otu.net <- cutreeDynamic(otu.tre, minClusterSize = 20)
names(otu.net) <- colnames(lund.interesting)

## Put Taxonomies to the modules
net.tax <- cbind(tax[match(names(otu.net), row.names(tax)),], Module = otu.net)
net.tax <- subset(net.tax, Module != 0)

## Get methanogen containing modules
meth.tax <- c("Methanobacterium", "Methanocella", "Methanosarcina", "Methanosaeta")
meth.mods <- unique(net.tax[net.tax$Genus%in%meth.tax,]$Module)
meth.mods <- net.tax[net.tax$Module%in%meth.mods,]

## Check for enrichment of genera in methane networks compared to whole dataset

gen.enr <- function(sub.net, whole.net){
  Genus = NA
  pval = NA
  sub.len <- nrow(sub.net)
  whole.len <- nrow(whole.net)
  genera <- as.character(unique(sub.net$Genus))
  genera <- genera[!genera%in%"unclassified"]
  
  ## Loop through genera and find enrichment
  for(i in 1:length(genera)){
    genus <- genera[i]
    cat(genus)
    sub.gen.count <- length(sub.net$Genus[sub.net$Genus == genus])
    whole.gen.count <- length(whole.net$Genus[whole.net$Genus == genus])
    p.val <- 1- phyper(sub.gen.count, whole.gen.count, 
                    whole.len - whole.gen.count, sub.len)
    Genus = c(Genus, genus)
    pval = c(pval, p.val)
    
  }
  genera.stats <- data.frame(Genus = Genus, Pval = pval)
  genera.stats$Padj <- p.adjust(genera.stats$Pval, method = "BH")
  return(genera.stats)
}

gen.stats <- gen.enr(meth.mods, net.tax)
subset(gen.stats[order(gen.stats$Padj),], Padj <= 0.05)
