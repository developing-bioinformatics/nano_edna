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
  summarise (n = n()) %>%
  group_by(samname) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(freq))

#NMDS
cast.gen = reshape2::dcast(t1.lca.genus[,1:3], samname ~ genus)
cast.gen[is.na(cast.gen)] = 0
rownames(cast.gen) = cast.gen[,1]
grp = cast.gen %>% tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'year', 'qual'), sep='[_]') 
mMDS = metaMDS(cast.gen[,-1], try=1000, k=2, distance = 'bray')
plot(mMDS, display='sites', type = 'p')
with(grp, ordihull(mMDS, site, draw='polygon', label=TRUE))

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
mMDS.sam = metaMDS(cast.gen.sam[,-1], try=1000, k=2, distance = 'bray')
plot(mMDS.sam, display='sites', type = 'p')
with(grouping, ordiellipse(mMDS.sam, site, draw='polygon', label=TRUE))

adonis2(cast.gen.sam[,-1] ~ site, data = grouping)







cast.gen = reshape2::dcast(t1.lca.genus, samname ~ genus)
cast.gen[is.na(cast.gen)] = 0
rownames(cast.gen) = cast.gen[,1]
gen.dist <- vegdist(decostand(cast.gen[,-1], method='normalize'), diag=T, method='jaccard')
plot(hclust(gen.dist))
#gen.pv <- pvclust(t(decostand(cast.gen[,-1], method='log')), 
  #                method.hclust="ward", nboot=100, parallel=T)
gen.pv <- pvclust(t(cast.gen[,-1]), 
                  method.hclust="ward", nboot=100, parallel=T)

plot(gen.pv, print.num=F, print.pv = T)


#try with sample grouping

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

# 
# # cluster on last_common
# t1.lca.lc = t1.lca %>%
#   group_by(QueryID, samname) %>%
#   dplyr::slice(1) %>%
#   group_by(last_common, samname) %>%
#   summarise(count = n()) %>%
#   arrange(desc(count)) %>%
#   filter(!is.na(last_common))  
# 
# #%>%
# # filter(count > 5) 
# 
# cast.lc = reshape2::dcast(t1.lca.genus, samname ~ last_common)
# cast.lc[is.na(cast.lc)] = 0
# rownames(cast.lc) = cast.lc[,1]
# lc.dist <- vegdist(decostand(cast.lc[,-1], method='log'), diag=T, method='ward')
# plot(hclust(lc.dist))
# lc.pv <- pvclust(t(decostand(cast.lc[,-1], method='log')), 
#                   method.hclust="ward", nboot=100, parallel=T)
# plot(lc.pv, print.num=F, print.pv = T)


ggplot(data = t1.lca.genus %>% filter(freq>0.2)) +
  geom_col(aes(x=samname, y = freq, fill=genus)) +
  theme_minimal() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 


# tile plot heatmap


ggplot(data = t1.lca.genus %>% filter(freq>0.1)) +
  geom_tile(aes(x=samname, y = genus, fill=freq)) +
  scale_fill_viridis_b() +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 


#need to deal with control sample parsing
t1.lca.genus.sam = t1.lca.genus %>%
  mutate(samname = if_else(stringr::str_detect(samname, 'control'), stringr::str_replace(samname, 'control', 'control_S1'), samname)) %>%
  tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'year', 'qual'), sep='[_]') %>%
  tidyr::unite(sample, c('site', 'sam')) %>%
  group_by(sample,genus) %>%
  summarize(count = sum(n)) %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  arrange(desc(freq)) %>%
  filter(!is.na(genus))  

ggplot(data = t1.lca.genus.sam %>% filter(genus != '') %>% filter(freq>0.001)) +
  geom_tile(aes(x=sample, y = genus, fill=freq)) +
  scale_fill_viridis_b() +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 

## put tile plot alongside dendrogram of kmer clustering


cast.gen.sam = reshape2::dcast(t1.lca.genus.sam, sample ~ genus)
cast.gen.sam[is.na(cast.gen.sam)] = 0
rownames(cast.gen.sam) = cast.gen.sam[,1]
gen.dist <- vegdist(decostand(cast.gen.sam[,-1], method='normalize'), diag=T, method='jaccard')
plot(hclust(gen.dist))

#gen.pv <- pvclust(t(decostand(cast.gen.sam[,-1], method='log')), 
#                method.hclust="ward", nboot=100, parallel=T)
gen.pv <- pvclust(t(decostand(cast.gen.sam[,-1], method='normalize')), 
                  method.hclust="ward.D2", nboot=100, parallel=T)

plot(gen.pv, print.num=F, print.pv = T)


#other ideas

#merge samples + plot
t1.sample.genus = t1.lca %>%
  mutate(sample = str_trunc(samname, 10)) %>%
  group_by(QueryID, samname) %>%
  dplyr::slice(1) %>%
  group_by(genus, sample) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(genus))  %>%
  filter(count > 5) 



ggplot(data = t1.sample.genus %>% mutate(total = sum(count))) +
  geom_col(aes(x=sample, y = count, fill=genus), position=position_dodge2(width = 0.9, preserve = "single")) +
  scale_y_log10() +
  theme_minimal() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 


#merge samples + plot
t1.family = t1.lca %>%
  mutate(sample = str_trunc(samname, 13)) %>%
  group_by(QueryID, samname) %>%
  dplyr::slice(1) %>%
  group_by(family, samname) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(family))  %>%
  filter(count > 20) 

ggplot(data = t1.family %>% mutate(total = sum(count))) +
  geom_col(aes(x=samname, y = count, fill=family), position=position_dodge2(width = 0.9, preserve = "single")) +
  scale_y_log10() +
  theme_minimal() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 

cast.fam = reshape2::dcast(t1.family, samname ~ family)
cast.fam[is.na(cast.fam)] = 0
rownames(cast.fam) = cast.fam[,1]
fam.dist <- vegdist(decostand(cast.fam[,-1], method='log'), diag=T, method='jaccard')
plot(hclust(fam.dist))
fam.pv <- pvclust(t(decostand(cast.fam[,-1], method='log')), 
                  method.hclust="ward.D2", nboot=100, parallel=T)
plot(fam.pv, print.num=F, print.pv = T)



# add chart to end of dendrogram?
