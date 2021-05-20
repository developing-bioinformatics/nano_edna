# convert BLAST database to Kraken2 database
# run get_amplicon_db.R first
#source('get_amplicon_db.R')

buildKrakenDBfromBLASTdb <- function(blast.db = 'blast_db', kraken.db='kraken_db', threads=2){
  if(dir.exists(paste(kraken.db, "/taxonomy"))){
    
  } else{
    getKDBtaxonomy = paste('kraken2-build --download-taxonomy --db ', blast.db)
    #system(getKDBtaxonomy) # get taxonomy files
  }  
  
  fasta = list.files(blast.db, pattern='fasta$', full.names = T)
  print(fasta)
  
  for(f in fasta) {
    kdb_add = paste('kraken2-build --add-to-library', f, '--db', kraken.db)
    system(kdb_add)
  }
  kdb_build = paste('kraken2-build --build --threads', threads, '--db', kraken.db)
  system(kdb_build)
}


buildKrakenDBfromBLASTdb('blast_db', 'kraken_db', threads=nclus)
