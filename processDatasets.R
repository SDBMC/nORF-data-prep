#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Test for compatible argument. If not, return an error
if (length(args) == 1) {
  if (!(args[1] %in% c("openprotEnsembl", "openprotRefseq","sorfs", "sorfsMS", "all"))) {
    stop("Invalid argument. Choose one of (openprotEnsembl, openprotRefseq, sorfs, sorfsMS, all).", call.=FALSE)
  } 
} else {
  stop("Invalid argument. Choose one of (openprotEnsembl, openprotRefseq, sorfs, sorfsMS, all).", call.=FALSE)
}
source('processingFunctions.R')

library(tidyverse)
library(Biostrings)

#Function calls appropriate functions from 'processingFunctions.R' based on supplied argument
processDataset <- function(dataset) {
  if (dataset == "openprotEnsembl") {
    openprot <- processOpenProt(annotation = "ensembl")
    return(openprot)
  } else if (dataset == "openprotRefseq") {
    openprotRefseq <- processOpenProt(annotation = "refseq")
    return(openprotRefseq)
  } else if (dataset == "sorfs") {
    sorfs <- processSorfs()
    return(sorfs)
  } else if (dataset == "sorfsMS") {
    sorfsMS <- processSorfs(dataset = "sorfsMS")
    return(sorfsMS)
  } else if (dataset == "all") {
    openprot <- processOpenProt(annotation = "ensembl")
    sorfs <- processSorfs()
    combined <- combineNovelORFs(sorfs,openprot)
    return(combined)
  }}

novelORFtableOriginal <- processDataset(args[1])

novelORFtableNamed <- addIDs(novelORFtableOriginal)

#Convert to bed12 and write out
novelORFbed <- createBed12(novelORFtableNamed) 
write_tsv(novelORFbed, path = paste0(args[1],"_38.bed"), col_names = F)

#Convert to gtf and write out
novelORFgtf <- createGTF(novelORFtableNamed)
write.table(novelORFgtf, paste0(args[1],"_38.gtf"), col.names = F, row.names = F, sep = "\t", quote = F)


