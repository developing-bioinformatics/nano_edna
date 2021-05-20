library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(pvclust)
library(vegan)

## code to generate ordination figure

t1.lca.genus.sam = data.table::fread('t1.lca.split')
cast.gen.sam = reshape2::dcast(t1.lca.genus.sam[,1:3], sample ~ genus)
cast.gen.sam[is.na(cast.gen.sam)] = 0
rownames(cast.gen.sam) = cast.gen.sam[,1]
grouping = cast.gen.sam %>% tidyr::separate(sample, c('site', 'subsam')) %>% select(site, subsam)
mMDS.sam = metaMDS(cast.gen.sam[,-1], try=10000, k=3, distance = 'bray')

#Figure 3
png(filename = 'figures/ordiplot.png')
plot(mMDS.sam, display='sites', type = 't')
with(grouping, ordiellipse(mMDS.sam, site, draw='polygon', label=TRUE))
dev.off()

pdf('figures/ordiplot.pdf')
plot(mMDS.sam, display='sites', type = 't')
with(grouping, ordiellipse(mMDS.sam, site, draw='polygon', label=TRUE))
dev.off()

#stats
permutest = adonis2(wisconsin(cast.gen.sam[,-1]) ~ site, data = grouping, permutations=100000, method = 'bray')


#Figure 2: Clustering samples and heatmap









