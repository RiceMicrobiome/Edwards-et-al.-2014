field.map <- read.table("~/RMB/Publication/Data/FieldExp/field_map.txt", header = T, row.names = 1)
gh.map <- read.table("~/RMB/Publication/Data/GreenhouseExp/gh_map.txt", header = T, row.names = 1)

counts <- read.table("~/RMB/Publication/Data/whole_otu_table.txt", header = T, row.names = 1)

field.counts <- counts[,match(row.names(field.map), colnames(counts))]
field.counts <- field.counts[rowSums(field.counts) > 0,]
field.samp.counts <- cbind(field.map, SampleCounts = colSums(field.counts))
write.table(field.samp.counts, file = "~/RMB/Publication/Data/FieldExp/sample_counts.txt", quote = F, sep = "\t")

gh.counts <- counts[,match(row.names(gh.map), colnames(counts))]
gh.counts <- gh.counts[rowSums(gh.counts) > 0,]
gh.samp.counts <- cbind(gh.map, colSums(gh.counts))
write.table(gh.samp.counts, file = "~/RMB/Publication/Data/GreenhouseExp/sample_counts.txt", quote = F, sep = "\t")
