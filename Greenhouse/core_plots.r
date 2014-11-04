gh.da <- read.table("~/RMB/Publication/Data/GreenhouseExp/DA_comp_site_OTUs.txt", header = T, row.names = 1)
field.da <- read.table("~/RMB/Publication/Data/FieldExp/enriched_otus.txt", header = T, row.names = 1)
tax <- read.table("~/RMB/Publication/Data/FieldExp/field_tax.txt", header = T, row.names = 1)
## Get core greenhouse first
load("~/RMB/Publication/Data/GreenhouseExp/glm.gh.rda")
gh.da.e <- subset(gh.comp.site.glm, color == "E")
gh.e.arb <- as.character(subset(gh.comp.site.glm, color == "E" & Site == "Arbuckle" & padj < 0.01)$OTU)
gh.e.dav <- as.character(subset(gh.comp.site.glm, color == "E" & Site == "Davis" & padj < 0.01)$OTU)
gh.e.sac <- as.character(subset(gh.comp.site.glm, color == "E" & Site == "Sacramento" & padj < 0.01)$OTU)
gh.core <- intersect(gh.e.arb, intersect(gh.e.dav, gh.e.sac))

field.da.e <- subset(field.da)
b1.e <- unique(as.character(subset(field.da, Site == "BB P1" & padj < 0.01)$OTU))
b4.e <- unique(as.character(subset(field.da, Site == "BB P4" & padj < 0.01)$OTU))
d18.e <- unique(as.character(subset(field.da, Site == "Ditaler 18" & padj < 0.01)$OTU))
d19.e <- unique(as.character(subset(field.da, Site == "Ditaler 19" & padj < 0.01)$OTU))
ds.e <- unique(as.character(subset(field.da, Site == "DS RR" & padj < 0.01)$OTU))
sch.e <- unique(as.character(subset(field.da, Site == "Scheidec" & padj < 0.01)$OTU))
sft.e <- unique(as.character(subset(field.da, Site == "SFT 20 A" & padj < 0.01)$OTU))
sp.e <- unique(as.character(subset(field.da, Site == "Spooner Airstrip" & padj < 0.01)$OTU))
field.all <- table(c(b1.e, b4.e, d18.e, d19.e, ds.e, sch.e, sft.e, sp.e))
field.core <- Reduce(intersect, list(b1.e, b4.e, d18.e, d19.e, ds.e, sch.e, sft.e, sp.e))

all.core <- intersect(field.core, gh.core)
all.core.tax <- tax[match(all.core, row.names(tax)),]
ggplot(all.core.tax, aes(x = Class, fill = Order)) +
  geom_bar() +
  coord_flip() +
  theme(text = element_text(size = 20))
