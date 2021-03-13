library(ShortRead)
library(DECIPHER)
library(kmer)
library(ape)
library(parallel)
library(parallelDist)
library(dendextend)
library(ggplot2)
library(dplyr)
library(ggdendro)

nclus = 24
taxonomizr_path = '/usr/share/data/taxonomizr/'
files_in = 'data'
res_out = 'targdb_out'
source('R/core.R')

#download_sra(dir.out = 'data')

files = list.files(files_in, full.names = TRUE)

fq = readDNAStringSet(files[1], format='fastq')
reads = fq

parKcount = function(x){
  return(kcount(x, k=5))
}

cl = makeCluster(nclus, type = 'FORK')
splits = clusterSplit(cl, reads)
kcount = parLapply(cl, splits, function(x){kcount(as.DNAbin(x), k=5)})
stopCluster(cl)
kcount = lapply(kcount, function(x) as.data.frame(x))
kc = dplyr::bind_rows(kcount)
kdist = parDist(as.matrix(kc), method = 'canberra', threads = nclus)

hc = hclust(kdist)

# cutree splits
cuts = cutree(hc, h=155)

hc$labels = cuts
# recolor tips by label ID
# dend = as.dendrogram(hc)
# color_vec = as.factor(labels(dend, "labels"))
# dend = dend %>% color_branches(col=color_vec) %>% set("branches_lwd", 0.5)
# ggd1 = as.ggdend(dend)
# ggd1

# now get consensus sequence within cuts
consensusReads = list()
for(i in 1:max(cuts)){
for(i in 1:100){
  if(length(reads[cuts==i])>1){
    align = AlignSeqs(reads[cuts==i], iterations=10, processors=4)
    consensusReads[[i]] = ConsensusSequence(align,  threshold=0.5)

  }
  
}
  

writeFasta(consensusReads, file = 'tmp_consensus.fasta')

amplicon='18S'
t1.coll = BLAST_pipeline2(
  files[1],
  blast_args =  '-max_target_seqs 10000',
  blast_db = paste('blast_db/', amplicon, '.fasta', sep=''),
  tax_db = taxonomizr_path,
  parallel = T,
  nclus = nclus,
  save.hits = TRUE,
  E.max = 1e-50,
  Perc.Ident.min = 0,
  blast.type = 'blastn'
)
