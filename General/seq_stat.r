library(ggplot2)
library(reshape2)

## Set the working directory
setwd("~/RICE/Publication/Data/")

## Load in data
counts <- read.table("whole_otu_table.txt", header = T, row.names = 1)
bad.otus <- read.table("bad_otus", header = F)$V1
map <- read.table("soil_lund.map", header = T, row.names = 1, sep = "\t")

## Remove bad OTUs from counts
counts.good <- counts[!row.names(counts)%in%bad.otus,]

## Break into different experiments
field.map <- subset(map, Cultivation != "Greenhouse")
gh.map <- subset(map, Cultivation == "Greenhouse")
field.counts <- counts.good[,match(row.names(field.map), colnames(counts.good))]
gh.counts <- counts.good[,match(row.names(gh.map), colnames(counts.good))]

## Get statistics between the two fields
field.total.otu <- nrow(field.counts[apply(field.counts, 1, sum) > 0,])
gh.total.otu <- nrow(gh.counts[apply(gh.counts, 1, sum) > 0,])

## Get OTUs with counts greater than 5
field.good <- field.counts[apply(field.counts, 1, sum) > 5,]
gh.good <- gh.counts[apply(gh.counts, 1, sum) > 5,]
dim(field.good)
dim(gh.good)
