library(rBLAST)
library(ShortRead)
library(DECIPHER)
library(kmer)
library(ape)
library(msa)
library(parallel)
cores=24
blastdb = "/usr/share/data/ncbi/nt/nt.fa"


library(ShortRead)
library(kmer)
library(ape)
library(parallel)
library(parallelDist)
library(dendextend)
library(ggplot2)

nclus = 24
dna = readFastq('blasthits', pattern='Swamp_S1_psbA3')
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

cltax = data.table::fread('targdb_out/t1.coll.35.csv', nThread=2)
#run lca on cltax
source('https://raw.githubusercontent.com/developing-bioinformatics/eDNA_BLAST/master/R/core.R')
source('R/core.R')
lca.tax = lca(cltax, parallel=F, nclus = nclus)
lca.labels = lca.tax %>%
  group_by(QueryID) %>%
  slice(1) 
hc$labels = lca.labels$last_common %>% tidyr::replace_na('unclassified')

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








otu = otu(as.DNAbin(reads), k = 5, threshold = 0.9, method = "centroid", nstart = 20)


cid = unique(otu)
seqs.clusters = list()
for(i in 1:length(cid)){
  seqs.clusters[[i]] = reads[which(otu==cid[i])]
}

align_cons = function(x){
  cores = 32
  if (length(x) > 1) {
    #cluster reads
    seqs <- OrientNucleotides(x)
    # perform the alignment
    aligned <- AlignSeqs(seqs, processors = cores)
    otu.reads = ConsensusSequence(aligned, threshold = 0.9, ambiguity=F)
    return(otu.reads)
  } else {
    return(x)
  }
}
cons = lapply(seqs.clusters, function(x) align_cons(x))
cons.app = cons[[1]]
for(i in 2:length(cons)){
  cons.app = append(cons.app, cons[[i]])
}
cons.app

#^^ BLAST these






#cluster reads
seqs <- OrientNucleotides(reads)
# perform the alignment
aligned <- AlignSeqs(seqs, processors=cores)
dend_stag_df <- StaggerAlignment(aligned, processors=cores)
distanceMat <- DistanceMatrix(dend_stag_df, processors=cores) #compute distance matrix
#distanceMat <- DistanceMatrix(seqs, type = 'dist', processors=cores) #compute distance matrix
#dist2 <- kdistance(as.DNAbin(seqs), k = 5) 
clustersId <- IdClusters(distanceMat, 
                         method="complete",
                         cutoff=0.5,
                         showPlot = TRUE,
                         myXStringSet = seqs,
                         processors=cores)
## get consensus
cid = unique(clustersId$cluster)
# group seqs by cluster ID
seqs.clusters = list()
for(i in 1:length(cid)){
  seqs.clusters[[i]] = seqs[which(clustersId$cluster==cid[i])]
}
cons = lapply(seqs.clusters, function(x) ConsensusSequence(x, threshold=0.5))
cons.app = cons[[1]]
for(i in 2:length(cons)){
  cons.app = append(cons.app, cons[[i]])
}
cons.app
# blast this ^^

## phylogeny
#library(seqinr)
#library(ape)
#library(msa)
#alignment <- msa(cons.app, "Muscle")
# seq.align <- msaConvert(alignment, type='seqinr::alignment')
# d <- dist.alignment(seq.align, "identity")
# tree <- nj(d)
# plot(tree)


# raxml tree
#write phylip file
#alignment.file = 'out.phy'
#write.phylip(alignment, file = alignment.file)
#system('module load raxml')
#system(paste('raxmlHPC -p 100 -s ', alignment.file,' -m GTRCAT -n tree', sep = ''))


# run blast linear

# run blast with num_threads

# R parallel blast


