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


krep.files = list.files('kraken_rep', full.names = T)
l <- lapply(krep.files, data.table::fread)
for(i in 1:length(l)){
  filename = krep.files[[i]]
  samname = tools::file_path_sans_ext(basename(filename)) 
  l[[i]]$samname = samname
}
krep <- bind_rows(l)

# kproc <- krep %>% 
#   filter(V4=='G') %>% 
#   #filter(V1>=1) %>% 
#   group_by(samname) %>%
#   mutate(samname = if_else(stringr::str_detect(samname, 'control'), 
#                            stringr::str_replace(samname, 'control', 'Control_S1'), samname)) %>%
#   mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
#   mutate(samname = str_replace(samname, "_A_B", "_A"))
# 
# 
# #NMDS
# cast.gen = reshape2::dcast(kproc, formula=samname ~ V6, value.var='V1')
# cast.gen[is.na(cast.gen)] = 0
# rownames(cast.gen) = cast.gen[,1]
# grp = cast.gen %>% tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'replicate', 'year', 'qual'), sep='[_]') 
# mMDS = metaMDS(cast.gen[,-1], try=1000, k=2, distance = 'bray')
# plot(mMDS, display='site', type = 'p')
# with(grp, ordihull(mMDS, site, draw='polygon', label=TRUE))

kproc.sample = krep %>%
  filter(V4=='G') %>% 
  mutate(samname = if_else(stringr::str_detect(samname, 'control'), stringr::str_replace(samname, 'control', 'Control_S1'), samname)) %>%
  mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
  mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
  tidyr::separate(samname, c('num', 'site', 'sam', 'amplicon', 'replicate', 'year', 'qual'), sep='[_]') %>%
  tidyr::unite(sample, c('site', 'sam')) %>%
  filter(V1>=0.05) %>%
  group_by(sample,V6) %>%
  summarise(readcount = sum(V2)) %>%
  mutate(freq = readcount / sum(readcount)) %>% 
  filter(freq>0.01)

cast.gen.sam = reshape2::dcast(kproc.sample, sample ~ V6)
cast.gen.sam[is.na(cast.gen.sam)] = 0
rownames(cast.gen.sam) = cast.gen.sam[,1]
grouping = cast.gen.sam %>% tidyr::separate(sample, c('site', 'sam')) %>% select(site, sam)
mMDS.sam = metaMDS(cast.gen.sam[,-1], try=10000, k=3, distance = 'bray')
#Figure 3
png(filename = 'figures/kraken_ordiplot.png')
plot(mMDS.sam, display='sites', type = 't')
with(grouping, ordiellipse(mMDS.sam, site, draw='polygon', label=TRUE))
dev.off()

pdf('figures/kraken_ordiplot.pdf')
plot(mMDS.sam, display='sites', type = 't')
with(grouping, ordiellipse(mMDS.sam, site, draw='polygon', label=TRUE))
dev.off()

#stats
permutest = adonis2(wisconsin(cast.gen.sam[,-1]) ~ site, data = grouping, permutations=100000, method = 'bray')



# t1.lca.genus.bysite = t1.lca.split %>%
#   filter(phylum=='Streptophyta') %>%
#   filter(Perc.Ident>80) %>%
#   filter(genus!='') %>%
#   group_by(QueryID, site, sam_id, amplicon, replicate) %>%
#   dplyr::slice(1) %>%
#   group_by(genus, site, sam_id) %>%
#   summarise (n = n())  %>% 
#   group_by(site, sam_id) %>%
#   mutate(freq = n / sum(n)) %>%
#   filter(freq>0.05) %>%
#   arrange(desc(freq)) %>%
#   tidyr::unite(site_sam, c("site", "sam_id")) %>%
#   arrange(desc(site_sam))


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

khits_files = list.files('kraken_classif/', pattern='fastq', full.names = TRUE)
khits_files = khits_files[file.size(khits_files)>2000] #keep only files > 2000bytes
khits_fnam = tools::file_path_sans_ext(basename(khits_files))
khits_fq = fastqq(khits_files, k=5, khits_fnam)

khits_fq_data = bind_cols(khits_fq@probeLabel, khits_fq@nReads) 
names(khits_fq_data) = c('samname', 'raw_reads')
khits_fq_data = khits_fq_data %>%
  mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
  mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
  mutate(samname = str_replace(samname, "control", "control_S1")) %>%
  tidyr::separate(samname, c("samnum", "site", "sam_id", "amplicon", "replicate", "year", "minq"), sep = "_") %>%
  tidyr::unite(site_sam, c("site", "sam_id"))

# hit_files = list.files('blasthits/', pattern='fastq', full.names = TRUE)
# hit_files = hit_files[file.size(hit_files)>2000] #keep only files > 2000bytes
# hit_fnam = tools::file_path_sans_ext(basename(hit_files))
# hit_fq = fastqq(hit_files, k=5, hit_fnam)
# 
# hit_fq_data = bind_cols(hit_fq@probeLabel, hit_fq@nReads) 
# names(hit_fq_data) = c('samname', 'raw_reads')
# hit_fq_data = hit_fq_data %>%
#   mutate(samname = str_replace(samname, "_2019", "_B_2019")) %>%
#   mutate(samname = str_replace(samname, "_A_B", "_A")) %>%
#   mutate(samname = str_replace(samname, "control", "control_S1")) %>%
#   tidyr::separate(samname, c("samnum", "site", "sam_id", "amplicon", "replicate", "year", "minq"), sep = "_") %>%
#   tidyr::unite(site_sam, c("site", "sam_id"))
#   
pal = brewer.pal(name='Paired', n = 12)

ktileplot = ggplot(data = kproc.sample ) +
  geom_tile(aes(x=sample, y = V6, fill=freq)) +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.y = element_text(angle=45),
        axis.text.x = element_text(angle=45, vjust =1, hjust = 1)) +
  scale_fill_steps(low=pal[1], high=pal[2], n.breaks = 10) +
  xlab('') + ylab('')

kraw_readcount = ggplot(data = fq_data) +
  geom_col(aes(x=site_sam, y = raw_reads, fill= amplicon)) +
  theme_linedraw() +
  theme(legend.position = 'right',
        axis.text.y = element_text(angle=45),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0,200000) +
  ylab('Raw Reads') + 
  scale_fill_brewer(palette = 'Paired')
kraw_readcount

kclassif_readcount = ggplot(data = khits_fq_data) +
  geom_col(aes(x=site_sam, y = raw_reads, fill= amplicon)) +
  theme_linedraw() +
  theme(legend.position = 'none',
        axis.text.y = element_text(angle=45),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(0,200000) +
  ylab('Classified Reads') + 
  scale_fill_brewer(palette = 'Paired')
kclassif_readcount

kfinal = kraw_readcount + kclassif_readcount + ktileplot + plot_layout(ncol = 1, heights = c(1, 1, 3))
kfinal


ggsave('figures/figure2B.png', kfinal, height = 9, width=4)


# Cluster by sequence identity then rearrange tile plot
#download_sra(dir.out = 'data')
# files = list.files('kraken_classif/', pattern='fastq', full.names = TRUE)
# 
# 
# files = files[file.size(files)>2000] #keep only files > 2000bytes
# fnam = tools::file_path_sans_ext(basename(files))
# fq = fastqq(files, k=5, fnam)
# 
# #replace fq@filenames
# labels = tools::file_path_sans_ext( basename(fq@filenames) )
# labels = stringr::str_replace(labels, 'classified_', '')
# labels = stringr::str_replace(labels, '_2019_minq7', '')
# fq@probeLabel = labels
# set.seed(1234)
# fq.kmer = fq@kmer
# colnames(fq.kmer) = labels
# result <- pvclust(fq.kmer, 
#                   method.hclust="average", nboot=100, parallel=T)
# plot(result, print.num=F, print.pv = c('au'))
# 
# ggd1 <- as.ggdend(as.dendrogram((result)) %>%
#                     set('branches_k_color', k = 4) %>%
#                     set('branches_lwd', 0.6) %>%
#                     set('labels_colors', k = 4) %>%
#                     set('labels_cex', 0.4))
# (g1 = ggplot(ggd1, horiz=T, label_offset=4) + 
#     theme(axis.text.y = element_text(angle=45)))
# 
# ggsave(g1, file = 'test.png', height=4, width=6)
