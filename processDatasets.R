#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Test for compatible argument. If not, return an error
if (length(args) == 1) {
  if (!(args[1] %in% c("openprot", "openprot2pep","sorfs", "sorfsMS", "all"))) {
    stop("Invalid argument. Choose one of (openprot, openprot2pep, sorfs, sorfsMS, all).", call.=FALSE)
  } 
} else {
  stop("Invalid argument. Choose one of (openprot, openprot2pep, sorfs, sorfsMS, all).", call.=FALSE)
}
source('processingFunctions.R')

library(tidyverse)
library(Biostrings)

#Function calls appropriate functions from 'processingFunctions.R' based on supplied argument
processDataset <- function(dataset) {
  if (dataset == "openprot") {
    openprot <- processOpenProt()
    return(openprot)
  } else if (dataset == "openprot2pep") {
    openprot2pep <- processOpenProt(dataset = "openprot2pep")
    return(openprot2pep)
  } else if (dataset == "sorfs") {
    sorfs <- processSorfs()
    return(sorfs)
  } else if (dataset == "sorfsMS") {
    sorfsMS <- processSorfs(dataset = "sorfsMS")
    return(sorfsMS)
  } else if (dataset == "all") {
    openprot <- processOpenProt()
    sorfs <- processSorfs()
    combined <- combineNovelORFs(sorfs,openprot)
    return(combined)
  }}

novelORFtable <- processDataset(args[1])

#Convert to bed12 and write out
novelORFbed <- createBed12(novelORFtable)
write_tsv(novelORFbed, path = paste0(args[1],"_38.bed"), col_names = F)

#Convert to gtf and write out
novelORFgtf <- createGTF(novelORFtable)
write.table(novelORFgtf, paste0(args[1],"_38.gtf"), col.names = F, row.names = F, sep = "\t", quote = F)







