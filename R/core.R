## Functions for core BLAST pipeline


## Download SRA data files function
#
#
#

download_sra = function(proj="PRJNA605442", retmax=60, ... ){
  require(rentrez)
  require(ShortRead)
  require(Biostrings)
  require(xml2)
  require(dplyr)
  r_search <- entrez_search(db="sra", term=paste(proj, "[GPRJ]", sep=''), retmax=60)
  id_fetch <- entrez_fetch(db="sra", id=r_search$ids, rettype = 'xml')
  doc <- read_xml(id_fetch)
  SRAFile <- doc %>% 
    xml_find_all("//RUN_SET") %>% 
    xml_find_all("//SRAFile") %>% 
    xml_attr('url') 
  
  # ammend to allow searching on other terms e.g., amplicon gene target, sample name (Swamp, Pond), ????
  
  get_one = SRAFile[grepl("SRR11043468", SRAFile)]
  download.file(SRAFile[1], 'test.fastq')
  dna = readFastq('.', pattern='test.fastq') ## On to BLAST analysis >>>
  
  
  
}

## Lowest Common Ancestor
# To DO: -- Fix error with no taxonomy or single match
# Add filtering options?
# Add output options (e.g., level='genus')

lca2 = function(x) {
  require(dplyr)
  require(multidplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  numcount = countshnames==1
  lastuni = tail(names(shortnames[numcount==T]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  newtax <- as.list(ifelse(countshnames==1,shortnames,NA))
  
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni) %>%
    mutate(superkingdom = newtax$superkingdom) %>%
    mutate(phylum = newtax$phylum) %>%
    mutate(class = newtax$class) %>%
    mutate(order = newtax$order) %>%
    mutate(family = newtax$family) %>%
    mutate(genus = newtax$genus) %>%
    mutate(species = newtax$species)
  return(ret)
}
