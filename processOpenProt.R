#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
twoPeptide = F
# test if there is at least one argument: if not, return an error
if (length(args) == 1) {
  twoPeptide = T
  if (args[1] != "2peptide") {
    stop("Invalid argument", call.=FALSE)
  }
}
library(tidyverse)
library(Biostrings)

#Read in BED, FASTA, and TSV files
bedFile <- read_tsv("openProtAllPredicted_38.bed", col_names = F, col_types = "ciiciciicicc")
tsvFile <- read_tsv("openProtAllPredicted.tsv", skip = 1, col_types = "ccciddicciicccciiddcccd")
proteinFasta <- readAAStringSet("openProtAllPredicted_38.fasta") 

#Get IDs for subsets of interest: 2peptide and any evidence
#Also filters out some altProts that are isoforms on other transcripts and therefore unlikely to be true nORFs

IDalt2peptide <- tsvFile %>% 
  filter(`protein type` == "AltProt" & `MS score` > 1) %>% 
  filter(!grepl('II', `protein accession (others)`)) %>% 
  select(`protein accession numbers`) %>% 
  distinct()
IDaltEvidence <- tsvFile %>% 
  filter(`protein type` == "AltProt" & (`MS score` != 0 | `TE score` != 0)) %>% 
  filter(!grepl('II', `protein accession (others)`)) %>% 
  select(`protein accession numbers`) %>% 
  distinct() 


#Dataframe with AA seq
dfp <- data.frame(str_split(names(proteinFasta), "\\|", simplify = T)[,1], paste(proteinFasta))
colnames(dfp) <- c("altIDs", "AAseq")

  
gtfFile <- tibble(
  seqname = tsvFile$chr,
  source = "OpenProt",
  feature = "altORF",
  start = tsvFile$`start genomic coordinates`,
  end = tsvFile$`stop genomic coordinates`,
  score = ".",
  strand = tsvFile$strand,
  frame = ".",
  altIDs = tsvFile$`protein accession numbers`) 

if (twoPeptide) {
  #Create bed file
  BEDalt2peptide <- bedFile[which(bedFile$X4 %in% IDalt2peptide$`protein accession numbers`),]
  write_tsv(BEDalt2peptide, path = "openProtAlt2peptide_38.bed", col_names = F)
  #Create gtfFile
  gtfFile <- gtfFile %>% 
    filter(altIDs %in% IDalt2peptide$`protein accession numbers`) %>%
    left_join(dfp, by = "altIDs") %>%
    distinct() %>%
    mutate(attribute = paste0('gene_id "', altIDs,'"; transcript_id "', altIDs, '.X"; AA_seq "', AAseq,'";' )) %>%
    select(-altIDs, -AAseq)
  write.table(gtfFile, file = "openProtAlt2peptide_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
  
} else {
  #Create bed file
  BEDaltEvidence <- bedFile[which(bedFile$X4 %in% IDaltEvidence$`protein accession numbers`),]
  write_tsv(BEDaltEvidence, path = "openProtAltEvidence_38.bed", col_names = F)
  #Create gtf file
  gtfFile <- gtfFile %>% 
    filter(altIDs %in% IDaltEvidence$`protein accession numbers`) %>%
    left_join(dfp, by = "altIDs") %>%
    distinct() %>%
    mutate(attribute = paste0('gene_id "', altIDs,'"; transcript_id "', altIDs, '.X"; AA_seq "', AAseq,'";' )) %>%
    select(-altIDs, -AAseq)
  write.table(gtfFile, file = "openProtAltEvidence_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
}
