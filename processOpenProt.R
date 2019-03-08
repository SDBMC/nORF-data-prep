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
tsvFile <- read.delim("openProtAllPredicted.tsv", header = T, skip = 1)
dnaFasta <- readDNAStringSet("openProtAllPredicted_38.dna.fasta") 
proteinFasta <- readAAStringSet("openProtAllPredicted_38.fasta") 

#Get IDs for subsets of interest: 2peptide and any evidence
IDalt2peptide <- unique(tsvFile[which(tsvFile$`protein.type` == "AltProt" & (tsvFile$`MS.score` > 1)),]$`protein.accession.numbers`)
IDaltEvidence <- unique(tsvFile[which(tsvFile$`protein.type` == "AltProt" & (tsvFile$`MS.score` != 0 | tsvFile$`TE.score` != 0)),]$`protein.accession.numbers`)

#Dataframes with DNA and AA seqs
df <- data.frame(str_split(names(dnaFasta), "\\|", simplify = T)[,1], paste(dnaFasta))
colnames(df) <- c("altIDs", "DNAseq")
dfp <- data.frame(str_split(names(proteinFasta), "\\|", simplify = T)[,1], paste(proteinFasta))
colnames(dfp) <- c("altIDs", "AAseq")

  
gtfFile <- tibble(
  seqname = tsvFile$chr,
  source = "OpenProt",
  feature = "altORF",
  start = tsvFile$start.genomic.coordinates,
  end = tsvFile$stop.genomic.coordinates,
  score = ".",
  strand = tsvFile$strand,
  frame = "0",
  altIDs = tsvFile$protein.accession.numbers)

if (twoPeptide) {
  #Create bed file
  BEDalt2peptide <- bedFile[which(bedFile$X4 %in% IDalt2peptide),]
  write_tsv(BEDalt2peptide, path = "openProtAlt2peptide_38.bed", col_names = F)
  #Create gtfFile
  gtfFile <- gtfFile %>% 
    filter(altIDs %in% IDalt2peptide) %>%
    left_join(df, by = "altIDs" ) %>%
    left_join(dfp, by = "altIDs") %>%
    distinct() %>%
    mutate(attribute = paste0('gene_id "', altIDs,'"; transcript_id "', altIDs, '.X" transcript_sequence "', DNAseq,'"; AA_seq "', AAseq,'"' )) %>%
    select(-altIDs, -DNAseq, -AAseq)
  write.table(gtfFile, file = "openProtAlt2peptide_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
  
} else {
  #Create bed file
  BEDaltEvidence <- bedFile[which(bedFile$X4 %in% IDaltEvidence),]
  write_tsv(BEDaltEvidence, path = "openProtAltEvidence_38.bed", col_names = F)
  #Create gtf file
  gtfFile <- gtfFile %>% 
    filter(altIDs %in% IDaltEvidence) %>%
    left_join(df, by = "altIDs" ) %>%
    left_join(dfp, by = "altIDs") %>%
    distinct() %>%
    mutate(attribute = paste0('gene_id "', altIDs,'"; transcript_id "', altIDs, '.X" transcript_sequence "', DNAseq,'"; AA_seq "', AAseq,'"' )) %>%
    select(-altIDs, -DNAseq, -AAseq)
  write.table(gtfFile, file = "openProtAltEvidence_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
}






