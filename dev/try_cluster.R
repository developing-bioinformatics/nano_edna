library(seqTools)
library(ggplot2)
library(ggdendro)
library(pvclust)
library(dendextend)

#source('https://raw.githubusercontent.com/developing-bioinformatics/eDNA_BLAST/master/R/core.R')
source('R/core.R')
#download_sra(dir.out = 'data')
files = list.files('blasthits', pattern='fastq', full.names = TRUE)
files = files[file.size(files)>2000] #keep only files > 1000bytes
fnam = tools::file_path_sans_ext(basename(files))
fq = fastqq(files, k=5, fnam)

#replace fq@filenames
labels = tools::file_path_sans_ext( basename(fq@filenames) )
labels = stringr::str_replace(labels, 'hits_', '')
labels = stringr::str_replace(labels, '_2019_minq7', '')
fq@probeLabel = labels
cb = cbDistMatrix(fq)
hc = hclust(as.dist(cb), method='single')
plot(hc)


# dhc <- as.dendrogram(hc)
# # Rectangular lines
# # ddata <- dendro_data(dhc, type = "rectangle") 
# # p <- ggplot(segment(ddata)) + 
# #   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
# #   coord_flip() + 
# #   scale_y_reverse(expand = c(0.2, 0)) + 
# #   theme_dendro()
# # p
# 
# 
# ggdendrogram(hc)



set.seed(1234)
fq.kmer = fq@kmer
colnames(fq.kmer) = labels
result <- pvclust(fq.kmer, 
                  method.hclust="average", nboot=100, parallel=T)
plot(result, print.num=F, print.pv = T)

dend %>%
  pvclust_show_signif(result) %>%
  plot()

ggd1 <- as.ggdend(as.dendrogram((result)) %>%
  set('branches_k_color', k = 4) %>%
  set('branches_lwd', 0.6) %>%
  set('labels_colors', k = 4) %>%
  set('labels_cex', 0.4))
(g1 = ggplot(ggd1, horiz=T, label_offset=4) + 
  theme(axis.text.y = element_text(angle=45)))

ggsave(g1, file = 'test.png', height=4, width=6)

?as.ggdend

