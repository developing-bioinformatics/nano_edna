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

download_sra(amplicon='ITS2', sample = 'control', dir.out = 'data')

files = list.files('data', pattern='.fastq', full.names = TRUE)

BLAST_pipeline(files[1], 
               blast_args =  NULL,
               blast_db = '/usr/share/data/ncbi/nt/nt.fa',
               tax_db = '/usr/share/data/taxonomizr/',
               parallel = TRUE,
               nclus = 48)

