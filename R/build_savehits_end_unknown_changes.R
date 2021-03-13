library(ggplot)
library(dplyr)
library(stringr)
library(data.table)

#t1 = data.table::fread('t1.csv', sep = ',')  #HUGE

#t1.lca = data.table::fread('t1.lca.csv', sep = ',')  #HUGE
t1.lca.files = list.files('out', pattern='t1.lca.coll', full.names = T)
l <- lapply(t1.lca.files, fread, sep=",")
t1.lca <- bind_rows(l)
colnames(t1.lca)[20] = 'samname'


## READ ^^ subfiles from 'out' folder 


#tlca = read.csv('tlca.csv')


#colnames(tlca)[3] = 'samname'


t1.lca.genus = t1.lca %>%
  group_by(QueryID, samname) %>%
  dplyr::slice(1) %>%
  group_by(genus, samname) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(genus))  

#%>%
# filter(count > 5) 

cast.gen = reshape2::dcast(t1.lca.genus, samname ~ genus)
cast.gen[is.na(cast.gen)] = 0
rownames(cast.gen) = cast.gen[,1]
gen.dist <- vegdist(decostand(cast.gen[,-1], method='log'), diag=T, method='jaccard')
plot(hclust(gen.dist))
gen.pv <- pvclust(t(decostand(cast.gen[,-1], method='log')), 
                  method.hclust="average", nboot=100, parallel=T)
plot(gen.pv, print.num=F, print.pv = T)


ggplot(data = t1.lca.genus %>% mutate(total = sum(count))) +
  geom_col(aes(x=samname, y = count, fill=genus), position='dodge') +
  scale_y_log10() +
  theme_minimal() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 



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
t1.sample.family = t1.lca %>%
  mutate(sample = str_trunc(samname, 13)) %>%
  group_by(QueryID, samname) %>%
  dplyr::slice(1) %>%
  group_by(family, sample) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(family))  %>%
  filter(count > 20) 

ggplot(data = t1.sample.family %>% mutate(total = sum(count))) +
  geom_col(aes(x=sample, y = count, fill=family), position=position_dodge2(width = 0.9, preserve = "single")) +
  scale_y_log10() +
  theme_minimal() +
  theme(legend.position = 'right',
        axis.text.x = element_text(angle=90)) 
# add chart to end of dendrogram?
