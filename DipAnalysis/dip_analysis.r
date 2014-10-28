library(ggplot2)
library(reshape2)
library(scales)

setwd("~/RMB/Publication/Data/Dip/")

counts <- read.table("dip_otu_table.txt", header = T, row.names = 1)
mito <- read.table("~/RMB/Reference/mito_otus.txt", header = F)$V1
plastid <- read.table("~/RMB/Reference/plastid_otus.txt", header = F)$V1

map <- read.table("dip.map", header = T, row.names = 1, sep = "\t")
map <- map[match(colnames(counts), row.names(map)),]

category <- rep("Microbial", nrow(counts))
category[row.names(counts)%in%c(mito, plastid)] <- "Organellar"


cat.counts <- aggregate(counts, data.frame(category), sum)
row.names(cat.counts) <- cat.counts[,1]
cat.counts <- t(cat.counts[,-1])
cat.ratio <- (cat.counts / rowSums(cat.counts)) * 100
cat.map <- map[match(row.names(cat.ratio), row.names(map)),]
cat.m <- melt(cbind(cat.map, cat.ratio))
cat.m.sum <- summarySE(cat.m, measurevar = "value", groupvars = c("SampleType", "variable"))
cat.m.sum$SampleType <- factor(cat.m.sum$SampleType, levels = c("0h Pre", "0h Post", "24h Post", "Soil"))

ggplot(subset(cat.m.sum, SampleType != "Soil"), aes(x = SampleType, y = value / 100, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = (value / 100) - (se / 100), ymax = (value / 100) + (se / 100)), width = 0.2, 
                position = position_dodge(0.9)) +
  scale_fill_manual(values = c("royalblue4", "orange")) +
  theme_classic() +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(x = "", y = "Percent of Reads", fill = "Read Type") +
  theme(text = element_text(size = 20))

t.test(subset(cat.m, variable == "Microbial" & SampleType == "0h Pre")$value, 
       subset(cat.m, variable == "Microbial" & SampleType == "0h Post")$value)

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
