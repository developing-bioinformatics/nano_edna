## Functions for core BLAST pipeline


## Download SRA data files function
#
download_sra = function(proj="PRJNA605442", 
                        dir.out = 'data', 
                        SRR = NULL, 
                        amplicon = NULL, 
                        sample = NULL, 
                        retmax=60, ... ){
  # always search on project ID
  # sometimes search on sample name (Control, Swamp, Pond)
  # sometimes search on amplicon
  require(rentrez)
  require(xml2)
  require(dplyr)
  r_search <- entrez_search(db="sra", term=paste(proj, "[GPRJ]", sep=''), retmax=60)
  id_fetch <- entrez_fetch(db="sra", id=r_search$ids, rettype = 'xml')
  doc <- read_xml(id_fetch)
  RUNSet <- doc %>% 
    xml_find_all("//RUN_SET") 
  RUNInfo <- RUNSet %>%
    xml_find_all("//RUN")
  POOL <- RUNSet %>%
    xml_find_all("//Member")
  SRAFile = RUNSet %>% 
    xml_find_all("//SRAFile")
  
  accession = RUNInfo %>%
    xml_attr('accession')
  url = SRAFile %>% 
    xml_attr('url')
  filename = SRAFile %>% 
    xml_attr('filename')
  sample_name = POOL %>%
    xml_attr('sample_name')
  
  #if not exist, create dir 'data' or proj
  
  if(!dir.exists(dir.out)){
    dir.create(dir.out)
  } 
  
  # ammend to allow searching on other terms e.g., amplicon gene target, sample name (Swamp, Pond), ????
  
  if(!is.null(SRR)) {
    # download file for SRR
    get_one = grepl(SRR, SRAFile) & grepl('fastq', SRAFile)
    targ_url = url[get_one]
    file = filename[get_one]
    download.file(targ_url, paste(dir.out, file, sep='/'), method = 'wget')
    #dna = readFastq('.', pattern='test.fastq') ## On to BLAST analysis >>>
  } else if (!is.null(amplicon) & !is.null(sample)){
    get_list = grepl(amplicon, sample_name) & grepl(sample, sample_name) & grepl('fastq', SRAFile)
    targ_url = url[get_list]
    file = filename[get_list]
    for(g in 1:length(targ_url)){
      download.file(targ_url[g], paste(dir.out, file[g], sep='/'))
    }
  } else if(!is.null(amplicon)) {
    get_list = grepl(amplicon, sample_name) & grepl('fastq', SRAFile)
    targ_url = url[get_list]
    file = filename[get_list]
    for(g in 1:length(targ_url)){
      download.file(targ_url[g], paste(dir.out, file[g], sep='/'), method = 'wget')
    }
  } else if(!is.null(sample)){
    get_list = grepl(sample, sample_name) & grepl('fastq', SRAFile)
    targ_url = url[get_list]
    file = filename[get_list]
    for(g in 1:length(targ_url)){
      download.file(targ_url[g], paste(dir.out, file[g], sep='/'), method = 'wget')
    } 
  } else {
    get_list = grepl('fastq', SRAFile)
    targ_url = url[get_list]
    file = filename[get_list]
    for(g in 1:length(targ_url)){
      download.file(targ_url[g], paste(dir.out, file[g], sep='/'), method = 'wget')
    } 
  }
}

## Lowest Common Ancestor
# To DO: -- Fix error with no taxonomy or single match (check)
# Add filtering options? (moved to filter function)
# Add output options (e.g., level='genus')

lca = function(results, parallel=FALSE, nclus=4) {
  require(dplyr)
  require(multidplyr)
  
  sub_lca = function(x) {
    nombre = vector()
    taxnames = c('superkingdom',
                 'phylum',
                 'order',
                 'family',
                 'genus',
                 'species')
    shortnames = apply(x[, taxnames], 2, unique)
    countshnames = sapply(shortnames, length)
    numcount = countshnames == 1
    lastuni = tail(names(shortnames[numcount == T]), n = 1)
    nombre = as.data.frame(x[1, which(colnames(x) == lastuni)])
    newtax <- as.list(ifelse(countshnames == 1, shortnames, NA))
    
   # if(length(nombre) == 0){ 
    if(all(is.na(newtax)) | length(nombre) == 0) {
      ret = x %>%
        mutate(last_common = NA) %>%
        mutate(level_lca = NA) %>%
        mutate(superkingdom = NA) %>%
        mutate(phylum = NA) %>%
        mutate(class = NA) %>%
        mutate(order = NA) %>%
        mutate(family = NA) %>%
        mutate(genus = NA) %>%
        mutate(species = NA)
    } else {
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
    }
    return(ret)
  }
  if(parallel == TRUE){
    cluster <- new_cluster(nclus)
    cluster_library(cluster, 'dplyr')
    cluster_copy(cluster, 'sub_lca')
    #split groups across multiple CPU cores
    prepare = results %>%
      group_by(QueryID) %>%
      partition(cluster) %>%  #split groups into cluster units
      do({sub_lca(.)}) %>%
      collect()
  } else {
    prepare = results %>%
      group_by(QueryID) %>%
      do({sub_lca(.)})
  }
  return(prepare)
}

BLAST_pipeline = function(fastq, 
                          blast_args =  NULL,
                          blast_db = NULL,
                          tax_db = NULL,
                          parallel = FALSE,
                          nclus = 4,
                          save.hits=FALSE
) { 
  require(rBLAST)
  require(ShortRead)
  require(Biostrings)
  require(dplyr)
  require(taxonomizr)
  
  # load taxonomizr databases
  tax_db_files = list.files(tax_db, full.names = TRUE)
  nodes = tax_db_files[grepl('nodes', tax_db_files)]
  names = tax_db_files[grepl('names', tax_db_files)]
  accession = tax_db_files[grepl('accessionTaxa', tax_db_files)]
  taxaNodes<-read.nodes.sql(nodes)
  taxaNames<-read.names.sql(names)
  
  # read fastq
  dna = readFastq(fastq)
  reads = sread(dna)
  qscores = quality(dna) 
  
  ## blast
  bl <- blast(db=blast_db)
  
  # optional parallel
  if (parallel == TRUE) {
    require(parallel)
    wpredict = function(x){
      return(predict(bl, x, , BLAST_args = blast_args))
      
    }
    clus = makeCluster(nclus, type ='FORK');
    splits = clusterSplit(clus, reads)
    p_cl = parLapply(clus, splits, wpredict)
    stopCluster(clus)
    
    cl = bind_rows(p_cl)
    
    
  } else {
    cl <- predict(bl, reads, BLAST_args = blast_args)
  }
  accid = as.character(cl$SubjectID) # accession IDs of BLAST hits
  #takes accession number and gets the taxonomic ID
  ids<-accessionToTaxa(accid, accession)
  #taxlist displays the taxonomic names from each ID #
  taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
  cltax=cbind(cl,taxlist)
  
  return(cltax)
  
}



BLAST_pipeline2 = function(fastq, 
                          blast_args =  NULL,
                          blast_db = NULL,
                          tax_db = NULL,
                          parallel = FALSE,
                          nclus = 4,
                          save.hits=FALSE, 
                          E.max=1, 
                          Perc.Ident.min = 0) { 
  require(rBLAST)
  require(ShortRead)
  require(Biostrings)
  require(dplyr)
  require(tidyr)
  require(taxonomizr)
  
  #for testing
  # blast_args =  NULL
  # blast_db = '/usr/share/data/ncbi/nt/nt.fa'
  # tax_db = '/usr/share/data/taxonomizr/'
  # parallel = TRUE
  # nclus = 48
  # save.hits=TRUE
  # #end testing setup
  
  if(save.hits==TRUE){
    if(dir.exists('blasthits')){} else {dir.create('blasthits')}
  }
  # load taxonomizr databases
  tax_db_files = list.files(tax_db, full.names = TRUE)
  nodes = tax_db_files[grepl('nodes', tax_db_files)]
  names = tax_db_files[grepl('names', tax_db_files)]
  accession = tax_db_files[grepl('accessionTaxa', tax_db_files)]
  taxaNodes<-read.nodes.sql(nodes)
  taxaNames<-read.names.sql(names)
  
  # read fastq
  filename = basename(fastq)
  dna = readFastq(fastq)
  
  ## blast
  bl <- blast(db=blast_db)
  
  # optional parallel
  if (parallel == TRUE) {
    require(parallel)
    wpredict = function(x){
      pr = predict(bl, x@sread, BLAST_args = blast_args)
      # get query IDs
      hits = pr$QueryID %>% as.data.frame() 
      colnames(hits) = 'QueryID'
      hits =  hits %>% separate(QueryID, c("query", "id"), "_")
      pr$QueryID = hits[,2]
      
      accid = as.character(pr$SubjectID) # accession IDs of BLAST hits
      #takes accession number and gets the taxonomic ID
      ids<-accessionToTaxa(accid, accession)
      #taxlist displays the taxonomic names from each ID #
      taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
      cltax=cbind(pr,taxlist)
      
      cltax = BLAST_filter(cltax, E.max = E.max, Perc.Ident.min = Perc.Ident.min)
      return(cltax)
    }
    clus = makeCluster(nclus, type ='FORK');
    splits = clusterSplit(clus, dna)
    p_cl = parLapply(clus, splits, wpredict)
    stopCluster(clus)
    
    #return(list(p_cl, splits))
    # do hit collection
    splits_keep = list()
    for(i in 1:length(splits)){
      splits_keep[[i]] = splits[[i]][unique(as.numeric(p_cl[[i]]$QueryID))] # keep only reads that have at least one blast hit
    }
    keep = splits_keep[[1]]
    for(z in 2:length(splits_keep)){
      keep = append(keep, splits_keep[[z]])
    }
    outfile = paste('blasthits/hits_', filename, sep='')
    if(file.exists(outfile)) { file.remove(outfile) }
    writeFastq(keep, file=outfile)
    
    
    #return(p_cl)
    
    #add split_ID to Query_ID for lca to work
    for(k in 1:length(p_cl)){
      if(nrow(p_cl[[k]]) == 0){ next } else {
        p_cl[[k]]$QueryID = paste(k, p_cl[[k]]$QueryID, sep='_')
      }
    }
    isnot = lapply(p_cl, nrow) > 0
    p_cl = p_cl[isnot]
    cl = bind_rows(p_cl)
    return(cl)
    
    
  } else {
    cl <- predict(bl, reads, BLAST_args = blast_args)
    hits = cl$QueryID %>% as.data.frame() 
    colnames(hits) = 'QueryID'
    hits =  hits %>% separate(QueryID, c("query", "id"), "_")
    cl$QueryID = hits[,2]
    
    accid = as.character(cl$SubjectID) # accession IDs of BLAST hits
    #takes accession number and gets the taxonomic ID
    ids<-accessionToTaxa(accid, accession)
    #taxlist displays the taxonomic names from each ID #
    taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
    cltax=cbind(cl,taxlist)
    cltax = BLAST_filter(cltax, E.max = E.max, Perc.Ident.min = Perc.Ident.min)
    
    #write hits file
    keep = dna[unique(as.numeric(cltax$QueryID))]
   
    outfile = paste('blasthits/hits_', filename, sep='')
    if(file.exists(outfile)) { file.remove(outfile) }
    writeFastq(keep, file=outfile)
    
    return(cltax)
  }

}


# Filter function
#f1 = BLAST_filter(t1, E.max= 1e-50)
BLAST_filter <- function(results = NULL, E.max=1, Perc.Ident.min = 0, ... ){
  require(dplyr)
  ret = results %>%
    filter(E <= E.max) %>%
    filter(Perc.Ident >= Perc.Ident.min)
  
  return(ret)
  
}

