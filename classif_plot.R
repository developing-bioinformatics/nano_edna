library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(pvclust)
library(vegan)
library(seqTools)
library(dendextend)
library(ggdendro)
library(RColorBrewer)

if(!file.exists('t1.lca.split')) {
  t1.lca.files = list.files('targdb_out/', pattern = 't1.coll', full.names = T)
  l <- lapply(t1.lca.files, fread, sep = ",")
  l2 = lapply(l, function(x) {
    x$QueryID = as.character(x$QueryID)
    return(x)
  }) # fix if QueryID types do not match
  t1.lca <- bind_rows(l2)
  rm(l)
  rm(l2)
  colnames(t1.lca)[21] = 'samname'
  
  
  t1.lca.split = t1.lca %>%
    mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
    mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
    mutate(samname = str_replace(samname, "control", "control_S1")) %>%
    tidyr::separate(
      samname,
      c(
        "samnum",
        "site",
        "sam_id",
        "amplicon",
        "replicate",
        "year",
        "minq"
      ),
      sep = "_"
    )
  
  #add write for this table to save a checkpoint
  fwrite(t1.lca.split, file = 't1.lca.split', nThread = 4)
  
} else {
  t1.lca.split = fread('t1.lca.split', nThread = 4)
}

t1.lca.genus.bysite = t1.lca.split %>%
  filter(phylum=='Streptophyta') %>%
  filter(Perc.Ident>80) %>%
  filter(genus!='') %>%
  group_by(QueryID, site, sam_id, amplicon, replicate) %>%
  dplyr::slice(1) %>%
  group_by(genus, site, sam_id, amplicon, replicate) %>% 
  mutate(count = n()) %>% slice(1) %>%
  group_by(site, sam_id, amplicon, replicate) %>% 
  mutate(total = sum(count)) %>%
  mutate(freq = count/total) %>%
  filter(freq>=0.05) %>% 
  group_by(genus, site, sam_id) %>%
  mutate(bysite_count = sum(count)) %>%  slice(1) %>% 
  tidyr::unite(site_sam, c("site", "sam_id")) %>%
  group_by(site_sam) %>%
  mutate(bysite_total = sum(bysite_count)) %>%
  mutate(bysite_freq = bysite_count / bysite_total) %>% 
  filter(bysite_freq>0.01) %>%
  arrange(desc(freq)) %>%
  arrange(desc(site_sam))


# get total read counts for each sample and calculate % classified
files = list.files('data', pattern='fastq', full.names = TRUE)
files = files[file.size(files)>2000] #keep only files > 2000bytes
fnam = tools::file_path_sans_ext(basename(files))
fq = fastqq(files, k=5, fnam)

fq_data = bind_cols(fq@probeLabel, fq@nReads) 
names(fq_data) = c('samname', 'raw_reads')
fq_data = fq_data %>%
  mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
  mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
  mutate(samname = str_replace(samname, "control", "control_S1")) %>%
  tidyr::separate(samname, c("samnum", "site", "sam_id", "amplicon", "replicate", "year", "minq"), sep = "_") %>%
  tidyr::unite(site_sam, c("site", "sam_id"))

hit_files = list.files('blasthits/', pattern='fastq', full.names = TRUE)
hit_files = hit_files[file.size(hit_files)>2000] #keep only files > 2000bytes
hit_fnam = tools::file_path_sans_ext(basename(hit_files))
hit_fq = fastqq(hit_files, k=5, hit_fnam)

hit_fq_data = bind_cols(hit_fq@probeLabel, hit_fq@nReads)
names(hit_fq_data) = c('samname', 'raw_reads')
hit_fq_data = hit_fq_data %>%
  mutate(samname = str_replace(samname, "hits_", "")) %>%
  mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
  mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
  mutate(samname = str_replace(samname, "control", "control_S1")) %>%
  tidyr::separate(samname, c("samnum", "site", "sam_id", "amplicon", "replicate", "year", "minq"), sep = "_") %>%
  tidyr::unite(site_sam, c("site", "sam_id"))
#   
pal = brewer.pal(name='Paired', n = 12)

tileplot = ggplot(data = t1.lca.genus.bysite) +
  geom_tile(aes(x=site_sam, y = genus, fill=freq)) +
  theme_linedraw() +
  theme(legend.position = 'none',
        axis.text.y = element_text(angle=45),
        axis.text.x = element_text(angle=45, vjust =1, hjust = 1)) + 
  scale_fill_steps(low=pal[1], high=pal[2], n.breaks = 10) +
  xlab('') + ylab('')
tileplot 

raw_readcount = ggplot(data = fq_data) +
  geom_col(aes(x=site_sam, y = raw_reads, fill= amplicon)) +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.y = element_text(angle=45),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, vjust=1, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer(palette = 'Paired') +
  ylim(0,200000) +
  ylab('Raw Reads')
raw_readcount


classif_readcount = ggplot(data = hit_fq_data) +
  geom_col(aes(x=site_sam, y = raw_reads, fill=amplicon)) +
  theme_linedraw() +
  theme(legend.position = 'none',
        axis.text.y = element_text(angle=45),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0,200000) +
  ylab('Classified Reads') + 
  scale_fill_brewer(palette = 'Paired')
classif_readcount

final = raw_readcount + classif_readcount + tileplot + plot_layout(ncol = 1, heights = c(1, 1, 3))

ggsave('figures/figure2A.png', final, height = 9, width=4)

source('R/kraken_savehits_end_targdb.R')

design <- "
  #11#
  2233
  4455
"

fig2_all = raw_readcount + 
  classif_readcount + 
  kclassif_readcount + 
  tileplot + 
  ktileplot + 
  plot_layout(design=design, heights=c(1,1,3))
fig2_all

ggsave('figures/figure2_all.png',fig2_all, height=14, width=7)



#ordiplot

cast.gen.sam = reshape2::dcast(t1.lca.genus.bysite %>% select(site_sam, genus, bysite_freq), site_sam ~ genus)
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
