#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(Guitar)
#Refseq entries


refseq <- makeTxDbFromUCSC(genome = "hg38", tablename = "refGene")
refseq <- renameSeqlevels(refseq , gsub("chr","", seqlevels(refseq)))
refseq <- keepStandardChromosomes(refseq )
source("classifynORFsConstants.R")

exonsByTranscript <- exonsBy(refseq, by = "tx", use.names = T)
intronsByTranscript <- intronsByTranscript(refseq, use.names = T)


#2.Read in novel ORFs list
novelORFs <- BED12toGRangesList("openprotRefseq_38.bed")
novelORFs <- renameSeqlevels(novelORFs, gsub("chr","", seqlevels(novelORFs)))
novelORFs <- keepStandardChromosomes(novelORFs)

transcriptTypes <- tibble(transcript_ID = exonsByTranscript@partitioning@NAMES, transcriptBiotype = NA) %>% 
  mutate(transcriptBiotype = ifelse(str_detect(transcript_ID, "NM_"), "protein_coding", transcriptBiotype)) %>% 
  mutate(transcriptBiotype = ifelse(str_detect(transcript_ID, "NR_"), "ncRNA", transcriptBiotype))


##### Set priority ranks for coding and non-coding annotations
#Priority rank for coding annotations
proteinCodingAnnotationRank <- tibble(ORFannotation = c("cds", "utr5-cds", "cds-utr3", "utr5", "utr3", "cds-intronic",
                                                        "utr5-intronic", "utr3-intronic", "cds-intergenic", "utr5-intergenic", "utr3-intergenic"),
                                      rank = 1:11)
#Priority rank for non-coding annotations
ncAnnotationRank <- tibble(transcriptBiotype = 'ncRNA', rank = 1)


proteinCodingTranscriptIDs <- transcriptTypes %>% 
  filter(transcriptBiotype == 'protein_coding')
proteinCodingExons <- exonsByTranscript[exonsByTranscript@partitioning@NAMES %in% proteinCodingTranscriptIDs$transcript_ID]
proteinCodingIntrons <- intronsByTranscript[intronsByTranscript@partitioning@NAMES %in% proteinCodingTranscriptIDs$transcript_ID]

nonCodingTranscriptIDs <- transcriptTypes %>%   
  filter(transcriptBiotype != 'protein_coding') 
nonCodingExons <- exonsByTranscript[exonsByTranscript@partitioning@NAMES %in% nonCodingTranscriptIDs$transcript_ID]
nonCodingIntrons <- intronsByTranscript[intronsByTranscript@partitioning@NAMES %in% nonCodingTranscriptIDs$transcript_ID]

#Separate nORFs into seven major classes and create tibble for annotating each novel ORF
#1. Within a protein coding exon
#2. Within a non coding exon
#3. Partial overlap of protein coding exon
#4. Partial overlap of non coding exon
#5. Within an intron of protein coding transcript
#6. Within an intron of non-coding transcript
#7. None of 1-3 (intergenic)
annotationClasses <- groupClasses()

#Classify in full detail
annotationMaster <- classify_nORFs(annotationTibble = annotationClasses, txdbInput = refseq)
annotationMasterFiltered <- annotationMaster %>% 
  mutate(codingRegion = ifelse(transcriptClass == "protein_coding" | transcriptClass == "intronic_codingTranscript" , T, F)) %>% 
  mutate(codingRegion = ifelse(transcriptClass == "intergenic", F, codingRegion)) %>% 
  dplyr::select(-transcriptClass) %>% 
  arrange(novelORF_ID)
write_tsv(annotationMasterFiltered, "nORF_refseqClassification.tsv")
