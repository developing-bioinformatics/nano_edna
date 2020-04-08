# test core functions
library(Biostrings)
library(taxonomizr)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringi)
library(forcats)
library(rBLAST)


source('R/core.R')

download_sra(amplicon='rbcLa', sample = 'Swamp', dir.out = 'data')

files = list.files('data', pattern='.fastq', full.names = TRUE)



t1 = BLAST_pipeline(files[8], 
               blast_args =  NULL,
               blast_db = '/usr/share/data/ncbi/nt/nt.fa',
               tax_db = '/usr/share/data/taxonomizr/',
               parallel = TRUE,
               nclus = 48)

# filter

f1 = BLAST_filter(t1, E.max= 1e-50)

f1.lca = lca(f1, parallel=T,  nclus = 24)


# plot
f1.lca %>% 
  group_by(QueryID) %>%
  dplyr::slice(1) %>%
  group_by(last_common) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

