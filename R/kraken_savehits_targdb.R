# test core functions
library(Biostrings)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringi)
library(forcats)
library(ggsci)

nclus = 32
files_in = 'data'
res_out = 'kraken_out'
kraken_reports = 'kraken_rep'
kreads_out = 'kraken_classif'
kdb = 'kraken_db'
source('R/core.R')

#download_sra(dir.out = 'data')

files = list.files(files_in, full.names = TRUE)

if(dir.exists(res_out)){} else {dir.create(res_out)}

for (i in 1:length(files)) {
  
  samname = tools::file_path_sans_ext(basename(files[[i]]))
  
  print(paste("Starting file:", samname, 'round', i, sep = ' '))
  print(date())
  p = proc.time()
  
  #files should be one or more fastq files
  
  #grep amplicon databases in ./blast_db
  if(grepl('control', samname)){
    amplicon = strsplit(samname, "_")[[1]][3]
  } else {
    amplicon = strsplit(samname, '_')[[1]][4]
  }
  if (grepl('psbA', amplicon)) {amplicon='psbA'}
  if (grepl('rbcL', amplicon)) {amplicon='rbcL'}
  
  print(amplicon)
  kreport = paste(kraken_reports, "/", samname, '.csv', sep = '')
  kout = paste(res_out, '/kraken.out.', samname, '.', amplicon, '.csv', sep = '')
  readsout = paste(kreads_out, '/classified.', samname, '.', amplicon, '.fastq', sep = '')
  #kraken2 --db $kdb  --threads 1 --use-names --report kreport.tab --fastq-input $targdata > kraken.out
  run_kraken = paste('kraken2 --db',
                     kdb,
                     '--threads',
                     nclus, 
                     '--use-names --report', 
                     kreport, 
                     '--classified-out',
                     readsout,
                     files[i],
                     ">", 
                     kout)
  system(run_kraken)
 
  xp = proc.time() - p
  print(paste("Last iteration took:", xp[[3]], sep = ' '))
}

#read e.g.,
# krep = data.table::fread('kraken_out/kraken.rep.9_Swamp_S2B_psbA3_2019_minq7.psbA.csv')
# krep %>% filter(V4=='G') %>% arrange(desc(V2))
