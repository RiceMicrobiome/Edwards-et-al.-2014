library(ggplot2)
library(reshape2)
library(VennDiagram)

## Load in the data
gh.otu <- read.table("~/RMB/Publication/Data/GreenhouseExp/gh_otu_table.txt", header = T, row.names = 1)
gh.map <- read.table("~/RMB/Publication/Data/GreenhouseExp/gh_map.txt", header = T, row.names = 1)
gh.da <- read.table("~/RMB/Publication/Data/GreenhouseExp/DA_OTUs.txt", header = T, row.names = 1)

load("~/RMB/TimeCourse/TC+GH/Data/tc_gh_data.rda")
otu.table <- tc.gh.data[[1]]
map <- tc.gh.data[[2]]
tax <- tc.gh.data[[3]]

tc.map <- subset(map, Cultivation == "TC")
tc.otu <- otu.table[,match(row.names(tc.map), colnames(otu.table))]
row.names(tc.otu) <- paste("Otu", row.names(tc.otu), sep = "")

## Get RS enriched only
RS.arb <- subset(gh.da, color == "RS" & Site == "Arbuckle")
RS.sac <- subset(gh.da, color == "RS" & Site == "Sacramento")
RS.dav <- subset(gh.da, color == "RS" & Site == "Davis")

## Get RP enriched only
RP.arb <- subset(gh.da, color == "RP" & Site == "Arbuckle")
RP.sac <- subset(gh.da, color == "RP" & Site == "Sacramento")
RP.dav <- subset(gh.da, color == "RP" & Site == "Davis")

## Get E enriched only
E.arb <- subset(gh.da, color == "E" & Site == "Arbuckle")
E.sac <- subset(gh.da, color == "E" & Site == "Sacramento")
E.dav <- subset(gh.da, color == "E" & Site == "Davis")

## RS depleted
RS.arb.down <- subset(gh.day)

### Venn Diagramms
venn.diagram(list( Endosphere = E.arb$OTU, Rhizoplane = RP.arb$OTU, Rhizosphere = RS.arb$OTU), 
             fill = c("#377EBA", "#4DAF4A", "#984EA3" ), cex = 3,
             filename = "~/RMB//Publication/GHFigures/DA/Site_Comp_Venn/arbuckle_comp.tiff")
venn.diagram(list( Endosphere = E.sac$OTU, Rhizoplane = RP.sac$OTU, Rhizosphere = RS.sac$OTU), 
             fill = c("#377EBA", "#4DAF4A", "#984EA3" ), cex = 3,
             filename = "~/RMB//Publication/GHFigures/DA/Site_Comp_Venn/sacramento_comp.tiff")
venn.diagram(list( Endosphere = E.dav$OTU, Rhizoplane = RP.dav$OTU, Rhizosphere = RS.dav$OTU), 
             fill = c("#377EBA", "#4DAF4A", "#984EA3" ), cex = 3,
             filename = "~/RMB//Publication/GHFigures/DA/Site_Comp_Venn/davis_comp.tiff")

arb.rp.only <- RP.arb$OTU[!RP.arb$OTU%in%c(as.character(RS.arb$OTU), as.character(E.arb$OTU))]
arb.e.only <- E.arb$OTU[!E.arb$OTU%in%c(as.character(RS.arb$OTU), as.character(RP.arb$OTU))]
sac.rp.only <- RP.sac$OTU[!RP.sac$OTU%in%c(as.character(RS.sac$OTU), as.character(E.sac$OTU))]
sac.e.only <- E.sac$OTU[!E.sac$OTU%in%c(as.character(RS.sac$OTU), as.character(RP.sac$OTU))]
dav.rp.only <- RP.dav$OTU[!RP.dav$OTU%in%c(as.character(RS.dav$OTU), as.character(E.dav$OTU))]
dav.e.only <- E.dav$OTU[!E.dav$OTU%in%c(as.character(RS.dav$OTU), as.character(RP.dav$OTU))]

core.rp <- intersect(arb.rp.only, intersect(dav.rp.only, sac.rp.only))
core.e <- intersect(E.arb$OTU, intersect(E.sac$OTU, E.dav$OTU))

### Put OTU table and map together
gh.whole <- cbind(gh.map, t(gh.otu))
dav.whole <- subset(gh.whole, Site == "Davis")
gh.whole$Compartment <- factor(gh.whole$Compartment, levels = c("Bulk Soil", "Rhizosphere","Rhizoplane", "Endosphere"))
dav.rp.only <- RP.dav$OTU[!RP.dav$OTU%in%E.dav$OTU]
pdf("~/Desktop/dav.rp.pdf")
for(x in 1:length(dav.rp.only)){
  otu = dav.rp.only[x]
  p <- ggplot(dav.whole, aes(x = Compartment, y = dav.whole[,colnames(dav.whole) == otu], fill = Compartment))
  p <- p + geom_boxplot()
  p <- p + scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA"))
  p <- p + labs(y = otu)
  print(p)
}
dev.off()

dav.e.only <- E.dav$OTU[!E.dav$OTU%in%c(as.character(RS.dav$OTU), as.character(RP.dav$OTU))]
pdf("~/Desktop/dav.e.pdf")
for(x in 1:length(dav.e.only)){
  otu = dav.e.only[x]
  p <- ggplot(dav.whole, aes(x = Compartment, y = dav.whole[,colnames(dav.whole) == otu], fill = Compartment))
  p <- p + geom_boxplot()
  p <- p + scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA"))
  p <- p + labs(y = otu)
  print(p)
}
dev.off()

ggplot(gh.whole, aes(x = Compartment, y = OtuNew.ReferenceOTU820, fill = Compartment)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  labs(y = otu) +
  facet_grid(.~Site)

tc.whole <- cbind(tc.map, t(tc.otu))
tc.whole$Compartment <- gsub("Bulk_Soil", "Bulk Soil", tc.whole$Compartment)
tc.whole$Compartment <- factor(tc.whole$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))

pdf("~/Desktop/e.core.pdf")
for (x in 1:length(core.e)){
  cat(otu)
  otu <- core.e[x]
  p <- ggplot(tc.whole, aes(x = "", y = tc.whole[,colnames(tc.whole) == otu], fill = Compartment))
  p <- p + geom_boxplot()
  p <- p + facet_grid(Compartment~Days, scales = "free_y")
  p <- p + scale_fill_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA"))
  p <- p + theme_bw()
  p <- p + theme(text = element_text(size = 20))
  p <- p + labs(y = otu)
  print(p)
}
dev.off()
  
tc.core.e.otu <- tc.otu[row.names(tc.otu)%in%core.e,]
tc.core.e.otu.quant <- tc.core.e.otu[apply(tc.core.e.otu, 1, sum) > 100,]

# check mean TMM
non.tc.core.e.otu <- tc.core.e.otu[!row.names(tc.core.e.otu)%in%row.names(tc.core.e.otu.quant),]
max(aggregate(t(non.tc.core.e.otu), data.frame(tc.map$Compartment), mean)[2,-1])
max(aggregate(t(tc.core.e.otu.quant), data.frame(tc.map$Compartment), mean)[2,-1])

tc.core.e.otu <- scale(t(tc.core.e.otu.quant), center = T)
tc.core.e.otu <- aggregate(tc.core.e.otu, data.frame(paste(tc.map$Compartment, tc.map$Days, sep = ".")), mean)
row.names(tc.core.e.otu) <- tc.core.e.otu[,1]
mdata <- colsplit(row.names(tc.core.e.otu), pattern = "\\.", names = c("Compartment", "Days"))
tc.core.e.otu <- cbind(mdata, tc.core.e.otu[,-1])
tc.core.e.otu$Days <- factor(tc.core.e.otu$Days)
tc.core.e.melt <- melt(tc.core.e.otu)
tc.core.e.melt$Compartment <- gsub("Bulk_Soil", "Bulk Soil", tc.core.e.melt$Compartment)
tc.core.e.melt$Compartment <- factor(tc.core.e.melt$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
tc.core.mean <- data.frame(tapply(tc.core.e.melt$value, paste(tc.core.e.melt$Compartment, tc.core.e.melt$Days, sep = "."), mean))
colnames(tc.core.mean) <- "value"
mdata <- colsplit(row.names(tc.core.mean), pattern = "\\.", names = c("Compartment", "Days"))
tc.core.mean <- cbind(mdata, tc.core.mean)
tc.core.mean$Days <- factor(tc.core.mean$Days)
tc.core.mean$Compartment <- factor(tc.core.mean$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))

tc.core.sd <- data.frame(tapply(tc.core.e.melt$value, paste(tc.core.e.melt$Compartment, tc.core.e.melt$Days, sep = "."), sd))
colnames(tc.core.sd) <- "value"
mdata <- colsplit(row.names(tc.core.sd), pattern = "\\.", names = c("Compartment", "Days"))
tc.core.sd <- cbind(mdata, tc.core.sd)
tc.core.sd$Days <- factor(tc.core.sd$Days)
tc.core.sd$Compartment <- factor(tc.core.sd$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Rhizoplane", "Endosphere"))

tc.core.mean13 <- subset(tc.core.mean, Days != 34 & Days != 55 & Days != 21)
tc.core.mean13$Days <- as.numeric(as.character(tc.core.mean13$Days))
tc.core.sd13 <- subset(tc.core.sd, Days != 34 & Days != 55 & Days != 21)
tc.core.sd13$Days <- as.numeric(as.character(tc.core.sd13$Days))


tc.summary <- summarySE(tc.core.e.melt, measurevar="value", groupvars=c("Days", "Compartment"))
tc.summary$Days <- as.numeric(as.character(tc.summary$Days))
pd <- position_dodge(0.2)
ggplot(subset(tc.summary, Days < 21), aes(group = Compartment, x = Days, y = value, color = Compartment)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.1, position = pd) +
  geom_line(size = 2, position = pd) +
  geom_point(position = pd) +
  theme_minimal() +
  scale_color_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  theme(axis.line = element_line(size = 1), axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 30), legend.key = element_blank()) +
  labs(y = "Scaled Abundance")
  
   

ggplot(subset(tc.core.e.melt, Days != 34 & Days != 55 & Days != 21), aes(x = as.numeric(as.character(Days)), y = value, group = variable, color = Compartment)) +
  geom_point(size = 1, position = position_jitter(width = 0.5)) +
  geom_line(aes(x = as.numeric(as.character(Days)), y = value, color = Compartment, group = Compartment), data = subset(tc.core.mean, Days != 34 & Days != 55 & Days != 21), inherit.aes = F, parse = F, size = 2) +
  #facet_grid(Compartment~.) +
  theme_bw() +
  scale_color_manual(values = c("#E41A1C", "#984EA3", "#4DAF4A", "#377EBA")) +
  labs(y = "Scaled Abundance") +
  theme(text = element_text(size = 20), legend.key = element_blank())
  

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}