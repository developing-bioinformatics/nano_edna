library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(pvclust)
library(vegan)

t1.lca.files = list.files('targdb_out/', pattern='t1.coll', full.names = T)
l <- lapply(t1.lca.files, fread, sep=",")
l2 = lapply(l, function(x) {x$QueryID = as.character(x$QueryID); return(x)}) # fix if QueryID types do not match
t1.lca <- bind_rows(l2)
colnames(t1.lca)[21] = 'samname'


t1.lca.genus = t1.lca %>%
  filter(phylum=='Streptophyta') %>%
  filter(genus!='') %>%
  group_by(QueryID, samname) %>%
  dplyr::slice(1) %>%
  group_by(genus, samname) %>%
  summarise (n = n())  %>% 
  group_by(samname) %>%
#  filter(sum(n)>10) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

# #NMDS
# cast.gen = reshape2::dcast(t1.lca.genus[,1:3], samname ~ genus)
# cast.gen[is.na(cast.gen)] = 0
# rownames(cast.gen) = cast.gen[,1]
# grp = cast.gen %>% tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'year', 'qual'), sep='[_]') 
# mMDS = metaMDS(cast.gen[,-1], try=1000, k=2, distance = 'bray')
# plot(mMDS, display='sites', type = 't')
# with(grp, ordihull(mMDS, site, draw='polygon', label=TRUE))

t1.lca.genus.sam = t1.lca.genus %>%
  mutate(samname = if_else(stringr::str_detect(samname, 'control'), stringr::str_replace(samname, 'control', 'Control_S1'), samname)) %>%
  tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'year', 'qual'), sep='[_]') %>%
  tidyr::unite(sample, c('site', 'sam')) %>%
  group_by(sample,genus) %>%
  summarize(count = sum(n)) %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  arrange(desc(freq)) %>%
  filter(!is.na(genus))  

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









