#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
massSpec = F
# test if there is at least one argument: if not, return an error
if (length(args) == 1) {
  massSpec = T
  if (args[1] != "massSpec") {
    stop("Invalid argument", call.=FALSE)
  }
}

library(tidyverse)

#Read in TXT file
if (massSpec) {
  rawFile <- read_tsv("sorfsDownload_MS.txt", col_types = 'cciiccciccccccc')
} else {
  rawFile <- read_tsv("sorfsDownload.txt", col_types = 'cciicccicccccc')
}
#Remove duplicates based on identical Chr, start, stop, strand, AND AAseq
rawFile <- rawFile %>%
  distinct(Chromosome, `Sorf start`, `Sorf end`, Strand, `AA-sequence` , .keep_all = TRUE) %>%
  filter(`Sorf start` >= 0)

#Convert splice info from sorfs.org into usable bed12 format
#First deal with simple no splice cases
noSplice <- rawFile %>%
  filter(is.na(`Spliced start parts`)) %>%
  mutate(blockCount = "1",
         blockSizes = `Sorf end` - (`Sorf start` - 1),
         blockStarts = "0")
#Cases with splicing
spliced <- rawFile %>%
  filter(!is.na(`Spliced start parts`))

blockCounter <- function(sorf) {
  blockCount <- length(str_split(sorf[6], "_", simplify = T))
  return(blockCount)
}
blockSizer <- function(sorf) {
  starts <- as.numeric(str_split(sorf[6], "_", simplify = T)) - as.numeric(sorf[3]) - 1
  ends <- as.numeric(str_split(sorf[7], "_", simplify = T)) - as.numeric(sorf[3])
  blockSizes <- paste((ends - starts), collapse = ",")
  return(blockSizes)
}
blockStarter <- function(sorf) {
  starts <- as.numeric(str_split(sorf[6], "_", simplify = T)) - as.numeric(sorf[3])
  blockStarts <- paste(starts, collapse = ",")
  return(blockStarts)
}

spliced <- spliced %>% 
  mutate(blockCount = apply(spliced, 1, blockCounter),
         blockSizes = apply(spliced, 1, blockSizer),
         blockStarts = apply(spliced, 1, blockStarter))

#Join back together
sorfsSpliceFormatted <- rbind(noSplice, spliced)

#Deal with similar sorfs by selecting longest when having the same end site + splice count

sorfsPlus <- sorfsSpliceFormatted %>%
  filter(Strand == '1') %>% 
  arrange(Chromosome,`Sorf end`, `Sorf start`) %>%
  distinct(Chromosome,`Sorf end`, blockCount, .keep_all = T)

sorfsMinus <- sorfsSpliceFormatted %>% 
  filter(Strand == '-1') %>% 
  arrange(Chromosome, `Sorf start`, -`Sorf end`) %>%
  distinct(Chromosome,  `Sorf start`, blockCount, .keep_all = T)

sorfsJoined <- bind_rows(sorfsPlus, sorfsMinus) %>% 
  arrange(Chromosome, `Sorf start`, `Sorf end`)

#Convert to bed12 format
bed12 <- tibble(
  chrom = paste0("chr", sorfsJoined$Chromosome),
  chromStart = sorfsJoined$`Sorf start` - 1,
  chromStop = sorfsJoined$`Sorf end`,
  name = gsub(":", "_", sorfsJoined$`Sorf ID`),
  score = "0",
  strand = gsub("1", "+", gsub("-1", "-", sorfsJoined$Strand)),
  thickStart = sorfsJoined$`Sorf start` - 1,
  thickEnd = sorfsJoined$`Sorf end`,
  itemRgb = "0,0,0",
  blockCount = sorfsJoined$blockCount,
  blockSizes = sorfsJoined$blockSizes,
  blockStarts = sorfsJoined$blockStarts)

#Remove a set of ~100 incorrectly annotated splice sites
bed12 <- bed12 %>%
  filter(substr(blockStarts, 1, 1) == "0")

#Create the gtf file
gtfFile <- tibble(
  seqname = sorfsJoined$Chromosome,
  source = "sorfs.org",
  feature = "sORF",
  start = sorfsJoined$`Sorf start`,
  end = sorfsJoined$`Sorf end`,
  score = ".",
  strand = gsub("1", "+", gsub("-1", "-", sorfsJoined$Strand)),
  frame = ".",
  attributes = paste0('gene_id "', sorfsJoined$`Sorf ID`,'"; transcript_id "', sorfsJoined$`Ensembl transcript ID`, 
                      '" transcript_sequence "', sorfsJoined$`Transcript sequence`,'"; AA_seq "', sorfsJoined$`AA-sequence`,
                      '"; start_codon "', sorfsJoined$`Start codon`, '"; sorf_length "', sorfsJoined$`Sorf length`,
                      '"; annotation "', sorfsJoined$Annotation, '"; biotype "', sorfsJoined$Biotype,'"' ))

#Write out bed and gtf files
if (massSpec) {
  write_tsv(bed12, path = "sorfsMS_38.bed", col_names = F)
  write.table(gtfFile, file = "sorfsMS_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
} else {
  write_tsv(bed12, path = "sorfs_38.bed", col_names = F)
  write.table(gtfFile, file = "sorfs_38.gtf", col.names = F, row.names = F, sep = "\t", quote = F)
}
