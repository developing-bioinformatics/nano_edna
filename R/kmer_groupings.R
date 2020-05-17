library(ShortRead)
library(kmer)
library(ape)
library(parallel)
library(parallelDist)
library(dendextend)
library(ggplot2)
library(dplyr)
library(ggdendro)

nclus = 24
dna = readFastq('blasthits', pattern='OHara_S1B_rbcL')
reads = sread(dna, id=id(dna))

parKcount = function(x){
  return(kcount(x, k=5))
}

cl = makeCluster(nclus, type = 'FORK')
splits = clusterSplit(cl, reads)
#kdist = kdistance(as.DNAbin(reads[1:100]), k = 5) # what if we parallelize kcount
kcount = parLapply(cl, splits, function(x){kcount(as.DNAbin(x), k=5)})
stopCluster(cl)
kcount = lapply(kcount, function(x) as.data.frame(x))
kc = dplyr::bind_rows(kcount)
kdist = parDist(as.matrix(kc), method = 'canberra', threads = nclus)

hc = hclust(kdist)
plot(hc)

cltax = data.table::fread('targdb_out/t1.coll.24.csv', nThread=2)
#run lca on cltax
source('https://raw.githubusercontent.com/developing-bioinformatics/eDNA_BLAST/master/R/core.R')
#source('R/core.R')
lca.tax = lca(cltax, parallel=T, nclus = nclus)
lca.labels = lca.tax %>%
  group_by(QueryID) %>%
  dplyr::slice(1) 
labels=  lca.tax %>% group_by(QueryID) %>% dplyr::slice(1) %>% ungroup() %>% select(family) 
labels = as.vector(labels$family) %>% tidyr::replace_na('unclassified')
hc$labels = labels
# recolor tips by label ID
dend = as.dendrogram(hc)
color_vec = as.factor(labels(dend, "label"))
dend = dend %>% color_branches(col=color_vec) %>% set("branches_lwd", 0.5)
ggd1 = as.ggdend(dend)
ggplot() +
  geom_segment(data = segment(ggd1), aes(x=x, y=y, xend=xend, yend=yend))+
  geom_segment(data = ggd1$segments %>%
                 filter(yend == 0) %>%
                 left_join(ggd1$labels, by = "x"), aes(x=x, y=y.x, xend=xend, yend=yend, color = label)) +
  coord_polar(theta="x") +
  scale_y_reverse(expand=c(0.2, 0)) +
  theme_dendro()
