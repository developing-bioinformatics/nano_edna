# script to get gene specific databases:
require(rentrez)
set_entrez_key("55549707f8ed138bcb423d98d739c7789408")
Sys.getenv("ENTREZ_KEY")


# to do list:
# add time reporting
# add flex query construction
# parallelize?


get_amplicon_db <- function(target, db.dir = 'blast_db', key = NULL) {
  #get fasta file with gene sequences from nuccore
  #write to blast database
  #works for plants only
  #no fuzzy matching of gene names? No handling of gene IDs. You get what you get.
  
  require(rentrez)
  require(rBLAST)

  by = 10000
  query = paste(target, '[All Fields] AND (plants[filter] AND ("0"[SLEN] : "100000"[SLEN]))',sep='') # allow more flexibility in searches
  genesearch <-
    entrez_search(db = "nuccore",
                  term = query,
                  use_history = TRUE, key=key)
  db_file = paste(db.dir, '/', target, ".fasta", sep='')
  if(dir.exists(db.dir)){} else {dir.create(db.dir)}
  get_gene = function(x) {
    recs <-
      entrez_fetch(
        db = "nuccore",
        web_history = genesearch$web_history,
        rettype = "fasta",
        retmax = by,
        retstart = x, 
        key=key
      )
    cat(recs, file = db_file, append = TRUE)
    cat(x + by, "sequences downloaded out of", genesearch$count, "\r")
  }
  sapp = lapply(seq(1, genesearch$count, by), get_gene)
  makeblastdb(db_file, dbtype = "nucl")
  return()
}

#get_amplicon_db('rbcL', key =Sys.getenv("ENTREZ_KEY"))
#get_amplicon_db('trnL', key =Sys.getenv("ENTREZ_KEY"))
#get_amplicon_db('psbA', key =Sys.getenv("ENTREZ_KEY"))
#get_amplicon_db('MATK', key =Sys.getenv("ENTREZ_KEY"))
#get_amplicon_db('ITS2', key =Sys.getenv("ENTREZ_KEY"))
#get_amplicon_db('18S', key =Sys.getenv("ENTREZ_KEY"))



