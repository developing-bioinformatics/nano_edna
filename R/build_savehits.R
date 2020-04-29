# test core functions
library(Biostrings)
library(taxonomizr)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringi)
library(forcats)
library(rBLAST)
library(ggsci)

nclus = 64

source('R/core.R')

#download_sra(dir.out = 'data')

files = list.files('data', full.names = TRUE)
t1.coll = list()
t1.lca.coll = list()
tlca.coll = list()

if(dir.exists('out')){} else {dir.create('out')}

for (i in 1:length(files)) {
  
  samname = tools::file_path_sans_ext(basename(files[[i]]))
  
  print(paste("Starting file:", samname, 'round', i, sep = ' '))
  print(date())
  p = proc.time()
  #files should be one or more fastq files
  t1.coll[[i]] = BLAST_pipeline2(
    files[i],
    blast_args =  NULL,
    blast_db = '/usr/share/data/ncbi/nt/nt.fa',
    tax_db = '/usr/share/data/taxonomizr/',
    parallel = TRUE,
    nclus = nclus,
    save.hits = TRUE,
    E.max = 1e-25,
    Perc.Ident.min = 0
  ) #filter on E value not Perc.Ident
  
  t1.coll[[i]] = cbind(t1.coll[[i]], rep(samname, nrow(t1.coll[[i]])))
  
  if(nrow(t1.coll[[i]]) == 0){ next }
  
  
  t1.lca.coll[[i]] = lca(t1.coll[[i]], parallel = T,  nclus = nclus)

  
  # plot
  tlca.coll[[i]] = t1.lca.coll[[i]] %>%
    group_by(QueryID) %>%
    dplyr::slice(1) %>%
    group_by(last_common) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  tlca.coll[[i]] = cbind(tlca.coll[[i]], rep(samname, nrow(tlca.coll[[i]])))
  
  write.csv(tlca.coll[[i]], paste('out/tlca.coll.', i, '.csv', sep = ''))
  write.csv(t1.lca.coll[[i]], paste('out/t1.lca.coll.', i, '.csv', sep = ''))
  write.csv(t1.coll[[i]], paste('out/t1.coll.', i, '.csv', sep = ''))
  
  xp = proc.time() - p
  print(paste("Last iteration took:", xp[[3]], sep = ' '))
}
#save.image('test.RData')
#load('test.RData')


t1 = bind_rows(t1.coll)

t1.lca = bind_rows(t1.lca.coll)

tlca = bind_rows(tlca.coll)

colnames(t1.lca.coll)[19] = 'samname'
colnames(tlca)[3] = 'samname'

write.csv(t1, 't1.csv')
write.csv(t1.lca, 't1.lca.csv')
write.csv(tlca, 'tlca.csv')

# ggplot(data = tlca) +
#   geom_col(aes(x=samname, y = count, col=last_common)) +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle=90))

