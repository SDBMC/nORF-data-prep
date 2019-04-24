#!/usr/bin/env Rscript

#Functions for processing sorfs and openProt datasets

processOpenProt <- function(annotation) {
  #Read in BED, FASTA, and TSV files based on ensembl or refseq
  
  if (annotation == "ensembl") {
    bedFile <- read_tsv("ensemblOpenProtAllPredicted_38.bed", 
                        col_names = c("chrom", "chromStart", "chromStop", "name", "score", 
                        "strand", "thickStart", "thickEnd", "itemRgb", "blockCount","blockSizes", "blockStarts"), 
                        col_types = "ciiciciicicc") %>% 
      #Remove a bugged entry
      filter(name != "IP_296985" | chrom != "chrY")
    
    tsvFile <- read_tsv("ensemblOpenProtAllPredicted.tsv", skip = 1, col_types = "ccciddicciicccciiddcccd")
    proteinFasta <- readAAStringSet("ensemblOpenProtAllPredicted_38.fasta") 
  } else if (annotation == "refseq") {
    
    bedFile <- read_tsv("refseqOpenProtAllPredicted_38.bed", 
                        col_names = c("chrom", "chromStart", "chromStop", "name", "score", 
                        "strand", "thickStart", "thickEnd", "itemRgb", "blockCount","blockSizes", "blockStarts"), 
                        col_types = "ciiciciicicc") 
    
    tsvFile <- read_tsv("refseqOpenProtAllPredicted.tsv", skip = 1, col_types = "ccciddicciicccciiddcccd")
    proteinFasta <- readAAStringSet("refseqOpenProtAllPredicted_38.fasta")     
      
    
  }
  

  
  #Get IDs for subsets of interest
  #Also filters out some altProts that are isoforms on other transcripts and therefore unlikely to be true nORFs
  
  IDaltEvidence <- tsvFile %>% 
    filter(`protein type` == "AltProt" & (`MS score` != 0 | `TE score` != 0)) %>% 
    filter(!grepl('II', `protein accession (others)`)) %>% 
    dplyr::select(`protein accession numbers`) %>% 
    distinct()
  
  #Dataframe with AA seq
  dfp <- data.frame(str_split(names(proteinFasta), "\\|", simplify = T)[,1], paste(proteinFasta))
  colnames(dfp) <- c("name", "aaSeq")

  #Create table formatted consistently with sorfs (effectively bed12 format + extra info)
  openProtTable <- bedFile %>% 
    filter(name %in% IDaltEvidence$`protein accession numbers`) %>% 
    mutate(source = "openprot.org") %>% 
    left_join(dfp, by = "name") %>%
    mutate(startCodon = "ATG")
  openProtTable <- openProtTable %>% 
    mutate(length = str_length(aaSeq))
  
  return(openProtTable)
}



processSorfs <- function(dataset = "sorfs") {
  #Read in TXT file
  if (dataset == "sorfsMS") {
    rawFile <- read_tsv("sorfsDownload_MS.txt", col_types = 'cciiccciccccccc')
  } else {
    rawFile <- read_tsv("sorfsDownload.txt", col_types = 'cciicccicccccc')
  }
  #Remove duplicates based on identical Chr, start, stop, strand, AND AAseq
  rawFile <- rawFile %>%
    distinct(Chromosome, `Sorf start`, `Sorf end`, Strand, `AA-sequence` , .keep_all = TRUE) %>%
    filter(`Sorf start` >= 0) %>% 
    mutate(`AA-sequence` = str_remove(`AA-sequence`, "\\*"))
  
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
    arrange(Chromosome, `Sorf start`, `Sorf end`) %>% 
    #Remove a set of ~100 incorrectly annotated splice sites
    filter(substr(blockStarts, 1, 1) == "0")
  
  #Create table formatted consistently with openProt (effectively bed12 format + extra info)
  sorfsTable <- tibble(
    chrom = paste0("chr", sorfsJoined$Chromosome),
    chromStart = sorfsJoined$`Sorf start` - 1,
    chromStop = sorfsJoined$`Sorf end`,
    name = gsub(":", "_", sorfsJoined$`Sorf ID`),
    score = 0,
    strand = gsub("1", "+", gsub("-1", "-", sorfsJoined$Strand)),
    thickStart = sorfsJoined$`Sorf start` - 1,
    thickEnd = sorfsJoined$`Sorf end`,
    itemRgb = "0,0,0",
    blockCount = as.integer(sorfsJoined$blockCount),
    blockSizes = sorfsJoined$blockSizes,
    blockStarts = sorfsJoined$blockStarts,
    source = "sorfs.org",
    aaSeq = sorfsJoined$`AA-sequence`,
    startCodon = sorfsJoined$`Start codon`,
    length = sorfsJoined$`Sorf length`)
  return(sorfsTable)
}

combineNovelORFs <- function(sorfs,openprot) {
  #Bind and remove duplcates between datasets
  merged <- bind_rows(sorfs,openprot) %>% 
    arrange(chrom, chromStart, chromStop, source) %>% 
    distinct(chrom, chromStart, chromStop, aaSeq, .keep_all = T)
  
  #Deal with similar orfs by selecting longest when having the same end site + splicing
  plusStrand <- merged %>%
    filter(strand == '+') %>% 
    arrange(chrom, chromStop, chromStart, source) %>%
    distinct(chrom, chromStop, blockCount, .keep_all = T)
  
  minusStrand <- merged %>% 
    filter(strand == '-') %>% 
    arrange(chrom, chromStart, -chromStop) %>%
    distinct(chrom, chromStart, blockCount, .keep_all = T)
  
  strandsJoined <- bind_rows(plusStrand, minusStrand) %>% 
    arrange(chrom, chromStart, chromStop)
  return(strandsJoined)
}

addIDs <- function(novelORFtable) {
  idKey <- read_csv("annotation_w_coordinates.csv", col_types = 'cccccccccccc')  %>% 
    mutate(mergeKey = str_c(V2,V3, c.original)) %>% 
    dplyr::select(mergeKey, V6)
  novelORFtableMerge <- novelORFtable %>% 
    mutate(mergeKey = str_c((chromStart +1),chromStop, name)) %>% 
    left_join(idKey, by = "mergeKey") %>% 
    mutate(name = V6) %>% 
    dplyr::select(-V6, mergeKey) %>% 
    filter(!(is.na(name)))
  return(novelORFtableMerge)
}

createBed12 <- function(novelORFtable) {
  return(novelORFtable[,1:12])
}

createGTF <- function(novelORFtable) {
  gtfFile <- tibble(
    seqname = gsub("chr", "", novelORFtable$chrom),
    source = novelORFtable$source,
    feature = "gene",
    start = novelORFtable$chromStart + 1,
    end = novelORFtable$chromStop,
    score = ".",
    strand = novelORFtable$strand,
    frame = ".",
    attributes = paste0('gene_id "', novelORFtable$name,'"; AA_seq "', novelORFtable$aaSeq,
                        '"; start_codon "', novelORFtable$startCodon, '"; sorf_length "', novelORFtable$length,'";' ))
  return(gtfFile)
}
