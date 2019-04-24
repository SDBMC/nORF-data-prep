#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(Guitar)

source("classifynORFsConstants.R")

##### First load in 3 required data files: txdb, novelORFs and transcript biotype
#1.Txdb for humans hg38
human <- makeTxDbFromEnsembl(organism="Homo sapiens",
                             release=96,
                             circ_seqs=DEFAULT_CIRC_SEQS,
                             server="ensembldb.ensembl.org",
                             username="anonymous", password=NULL, port=0L,
                             tx_attrib=NULL)


#2.Read in novel ORFs list
novelORFs <- BED12toGRangesList("all_38.bed")
novelORFs <- renameSeqlevels(novelORFs, gsub("chr","", seqlevels(novelORFs)))


#3.Load in transcript biotypes file
#In process group IG and TR genes into protein coding and re-annotate a few misclassified transcripts
#Note that transcriptTypes.txt file was generated for the ensembl biomart page with the follwing parameters:
  #Ensembl Genes 96, Dataset: Human genes(GRCh38.p12), Attributes: Transcript stable ID + Transcript type, Unique results only
proteinCodingBiotypes <- c("protein_coding","IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene",
                           "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")
transcriptTypes <- read_tsv("transcriptTypes.txt", col_names = c("transcript_ID", "transcriptBiotype"), skip = 1, col_types = "cc") %>% 
  mutate(transcriptBiotype = ifelse(transcriptBiotype %in% proteinCodingBiotypes, 'protein_coding', transcriptBiotype)) 


##### Set priority ranks for coding and non-coding annotations
#Priority rank for coding annotations
proteinCodingAnnotationRank <- tibble(ORFannotation = c("cds", "utr5-cds", "cds-utr3", "utr5", "utr3", "cds-intronic",
                                                        "utr5-intronic", "utr3-intronic", "cds-intergenic", "utr5-intergenic", "utr3-intergenic", "None"),
                                      rank = 1:12)
#Priority rank for non-coding annotations
ncAnnotationRank <- tibble(
  transcriptBiotype = c('nonsense_mediated_decay', 'non_stop_decay', 'retained_intron',
                        'sense_intronic', 'sense_overlapping', '3prime_overlapping_ncrna', 'antisense', 'bidirectional_promoter_lncRNA',
                        'IG_V_pseudogene', 'IG_C_pseudogene','TR_V_pseudogene', 'rRNA_pseudogene',
                        'polymorphic_pseudogene', 'translated_unprocessed_pseudogene', 'transcribed_processed_pseudogene',
                        'transcribed_unprocessed_pseudogene', 'processed_pseudogene', 'unprocessed_pseudogene', 
                        'unitary_pseudogene', 'transcribed_unitary_pseudogene', 'pseudogene', 'lincRNA', 
                        'macro_lncRNA', 'snoRNA','Mt_rRNA','Mt_tRNA','miRNA', 'snRNA','misc_RNA','processed_transcript', 'TEC'), 
  rank = 1:31)

exonsByTranscript <- exonsBy(human, by = "tx", use.names = T)
intronsByTranscript <- intronsByTranscript(human, use.names = T)
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
annotationMaster <- classify_nORFs(annotationTibble = annotationClasses, txdbInput = human)

#Label as protein coding T/F and in-frame with protein coding T/F
inFrameIDs <- generateInFrameIDs(gffFile = "Homo_sapiens.GRCh38.96.gff3", bed6File = "all_38.6.bed", txdb = human)
annotationMasterFiltered <- annotationMaster %>% 
  mutate(inFrame = ifelse(novelORF_ID %in% inFrameIDs$norf_ID, T, F)) %>% 
  mutate(codingRegion = ifelse(transcriptClass == "protein_coding" | transcriptClass == "intronic_codingTranscript" |
                                        transcriptBiotype == "nonsense_mediated_decay" | transcriptBiotype == "retained_intron", T, F)) %>% 
  mutate(codingRegion = ifelse(transcriptClass == "intergenic", F, codingRegion)) %>% 
  dplyr::select(-transcriptClass) %>% 
  arrange(novelORF_ID)
write_tsv(annotationMasterFiltered, "nORFclassification.tsv")
