
##### Specify functions used to classify norfs in classifynORFsCode.R
generateTranscriptPairs <- function(transcriptFeature, class, annotationTibble) {
  
  #Create subsetted Granges objects
  nORF_List <- annotationTibble %>%
    filter(transcriptClass == class)
  novelORFs <- novelORFs[novelORFs@elementMetadata@listData$V4 %in% nORF_List$novelORF_ID]
  #Find any and full overlaps
  anyoverlap <- findOverlaps(novelORFs, transcriptFeature, type = "any")
  transcriptPairsAnnotation <- tibble(novelORF_ID = novelORFs@elementMetadata@listData$V4[anyoverlap@from],
                                      transcript_ID = transcriptFeature@partitioning@NAMES[anyoverlap@to],
                                      transcriptClass = class) %>% 
    left_join(transcriptTypes, by = "transcript_ID") %>% 
    mutate(matchCompare = str_c(novelORF_ID, transcript_ID))
  #Check which ones are also full overlaps
  fulloverlap <- findOverlaps(novelORFs, transcriptFeature, type = "within")
  transcriptFullAnnotation <- tibble(novelORF_ID = novelORFs@elementMetadata@listData$V4[fulloverlap@from],
                                     transcript_ID = transcriptFeature@partitioning@NAMES[fulloverlap@to],
                                     transcriptClass = class) %>% 
    mutate(matchCompare = str_c(novelORF_ID, transcript_ID))
  
  #Label any matches that are also full matches
  transcriptPairsAnnotation <- transcriptPairsAnnotation %>% 
    mutate(isFull = ifelse(transcriptPairsAnnotation$matchCompare %in% transcriptFullAnnotation$matchCompare, T, F)) %>% 
    dplyr::select(-matchCompare)
  
  return(transcriptPairsAnnotation)
}


prioritizeNonCoding <- function(annotationPairs, ranks) {
  #Label the priority rank of each annotation
  annotationPairs <- annotationPairs %>% 
    left_join(ranks, by = "transcriptBiotype")
  #Sort by full overlaps first and then highest priority rank, then take first top result for each nORF
  annotationPairs <- annotationPairs %>% 
    arrange(novelORF_ID, -isFull, rank, transcript_ID) %>% 
    distinct(novelORF_ID, .keep_all = T) %>% 
    dplyr::select(-isFull, -rank)
  
  return(annotationPairs)
}

prioritizeCodingIntrons <- function(annotationPairs) {
  #Sort by full overlaps, then take first top result for each nORF
  annotationPairs <- annotationPairs %>% 
    arrange(novelORF_ID, -isFull, transcript_ID) %>% 
    distinct(novelORF_ID, .keep_all = T) %>% 
    dplyr::select(-isFull) 
  
  return(annotationPairs)
}

prioritizeCoding <- function(annotationPairs, ranks) {
  #Label the priority rank of each annotation
  annotationPairs <- annotationPairs %>% 
    left_join(ranks, by = "ORFannotation")
  #Sort by full overlaps first and then highest priority rank, then take first top result for each nORF
  annotationPairs <- annotationPairs %>% 
    arrange(novelORF_ID, rank, transcript_ID) %>% 
    distinct(novelORF_ID, .keep_all = T) %>% 
    dplyr::select(-rank)
  
  return(annotationPairs)
}


classifyProteinCoding <- function(annotationPairs, txdb) {
  
  #Create features
  introns <- intronsByTranscript(txdb, use.names = T)
  cds <- cdsBy(txdb, "tx", use.names = T)
  utr5 <- fiveUTRsByTranscript(txdb, use.names = T)
  utr3 <- threeUTRsByTranscript(txdb, use.names = T)
  
  annotOverlap <- function(novelORFs, annot, type) {
    overlap <- findOverlaps(novelORFs, annot, type = type)
    pairsTibble <- tibble(novelORF_ID = novelORFs@elementMetadata@listData$V4[overlap@from],
                          transcript_ID = annot@partitioning@NAMES[overlap@to]) %>% 
      mutate(matchCompare = str_c(novelORF_ID, transcript_ID))
  }
  
  #Find partial and full overlaps with features
  CDSany <- annotOverlap(novelORFs = novelORFs, annot = cds, type = "any")
  CDSfull <- annotOverlap(novelORFs = novelORFs, annot = cds, type = "within")
  UTR5any <- annotOverlap(novelORFs = novelORFs, annot = utr5, type = "any")
  UTR5full <- annotOverlap(novelORFs = novelORFs, annot = utr5, type = "within")
  UTR3any <- annotOverlap(novelORFs = novelORFs, annot = utr3, type = "any")
  UTR3full <- annotOverlap(novelORFs = novelORFs, annot = utr3, type = "within")
  INTRONany <- annotOverlap(novelORFs = novelORFs, annot = introns, type = "any")
  
  #Get comparable column
  annotationPairs <- annotationPairs %>% 
    mutate(matchCompare = str_c(novelORF_ID, transcript_ID))
  
  #Add feature overlaps for each transcript
  annotationPairs <- annotationPairs %>% 
    mutate(CDSany = ifelse(annotationPairs$matchCompare %in% CDSany$matchCompare, T, F)) %>%
    mutate(CDSfull = ifelse(annotationPairs$matchCompare %in% CDSfull$matchCompare, T, F)) %>%
    mutate(UTR5any = ifelse(annotationPairs$matchCompare %in% UTR5any$matchCompare, T, F)) %>%
    mutate(UTR5full = ifelse(annotationPairs$matchCompare %in% UTR5full$matchCompare, T, F)) %>%
    mutate(UTR3any = ifelse(annotationPairs$matchCompare %in% UTR3any$matchCompare, T, F)) %>%
    mutate(UTR3full = ifelse(annotationPairs$matchCompare %in% UTR3full$matchCompare, T, F)) %>%
    mutate(INTRONany = ifelse(annotationPairs$matchCompare %in% INTRONany$matchCompare, T, F)) %>%
    mutate(ORFannotation = "None") %>% 
    dplyr::select(-matchCompare)
  
  #Use feature overlaps to classify each norf-transcript pair
  annotationPairs <- annotationPairs %>% 
    mutate(ORFannotation = ifelse(CDSfull, "cds", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(CDSany & UTR5any, "utr5-cds", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(CDSany & UTR3any, "cds-utr3", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR5full, "utr5", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR3full, "utr3", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(CDSany & INTRONany & !CDSfull, "cds-intronic", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR5any & INTRONany, "utr5-intronic", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR3any & INTRONany, "utr3-intronic", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(CDSany & !CDSfull & !UTR3any & !UTR5any & !INTRONany, "cds-intergenic", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR5any & !UTR5full & !CDSany & !INTRONany, "utr5-intergenic", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(UTR3any & !UTR3full & !CDSany & !INTRONany, "utr3-intergenic", ORFannotation)) %>% 
    dplyr::select(novelORF_ID, transcript_ID, transcriptClass, ORFannotation) %>% 
    mutate(transcriptBiotype = "protein_coding")

  return(annotationPairs)
}

groupClasses <- function() {
  
  
  #1. Within a protein coding exon
  codingOverlaps <-  findOverlaps(novelORFs, proteinCodingExons, type = "within")
  codingOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(codingOverlaps@from)]
  
  #2. Within a non coding exon
  nonCodingOverlaps <- findOverlaps(novelORFs, nonCodingExons, type = "within")
  nonCodingOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(nonCodingOverlaps@from)]
  
  #3. Partial overlap of protein coding exon
  partialCodingOverlaps <-  findOverlaps(novelORFs, proteinCodingExons, type = "any")
  partialCodingOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(partialCodingOverlaps@from)]
  
  #4. Partial overlap of non coding exon
  partialnonCodingOverlaps <- findOverlaps(novelORFs, nonCodingExons, type = "any")
  partialnonCodingOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(partialnonCodingOverlaps@from)]
  
  #5. Within an intron of protein coding transcript
  codingIntronsOverlaps <-  findOverlaps(novelORFs, proteinCodingIntrons, type = "any")
  codingIntronsOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(codingIntronsOverlaps@from)]
  
  #6. Within an intron of non-coding transcript
  nonCodingIntronsOverlaps <- findOverlaps(novelORFs, nonCodingIntrons, type = "any")
  nonCodingIntronsOverlapsIDs <- novelORFs@elementMetadata@listData$V4[unique(nonCodingIntronsOverlaps@from)]
  
  
  #Create tibble for annotating each novel ORF
  annotationMaster <- tibble(novelORF_ID = novelORFs@elementMetadata@listData$V4,
                             transcript_ID = NA,
                             transcriptClass = NA,
                             transcriptBiotype = NA,
                             ORFannotation = NA) %>% 
    #Defaults intergenic
    mutate(transcriptClass = "intergenic") %>% 
    #Overwrite if non-coding intronic
    mutate(transcriptClass = ifelse(novelORF_ID %in% nonCodingIntronsOverlapsIDs, "intronic_ncTranscript", transcriptClass)) %>% 
    #Overwrite if coding intronic
    mutate(transcriptClass = ifelse(novelORF_ID %in% codingIntronsOverlapsIDs, "intronic_codingTranscript", transcriptClass)) %>% 
    #Overwrite if partial non-coding transcript exonic
    mutate(transcriptClass = ifelse(novelORF_ID %in% partialnonCodingOverlapsIDs, "non_coding", transcriptClass)) %>% 
    #Overwrite if partial coding transcript exonic
    mutate(transcriptClass = ifelse(novelORF_ID %in% partialCodingOverlapsIDs, "protein_coding", transcriptClass)) %>% 
    #Overwrite if full non-coding transcript exonic
    mutate(transcriptClass = ifelse(novelORF_ID %in% nonCodingOverlapsIDs, "non_coding", transcriptClass)) %>% 
    #Overwrite if full coding transcript exonic
    mutate(transcriptClass = ifelse(novelORF_ID %in% codingOverlapsIDs, "protein_coding", transcriptClass))
  return(annotationMaster)
}

classify_nORFs <- function(annotationTibble, txdbInput) {
  annotationTibble2 <- annotationTibble %>% 
    mutate(ORFannotation = ifelse(transcriptClass == "intergenic", "intergenic", ORFannotation))
  
  ncIntronsPairs <- generateTranscriptPairs(transcriptFeature = nonCodingIntrons, class = "intronic_ncTranscript", annotationTibble = annotationTibble2)
  annotation_ncIntrons <- prioritizeNonCoding(annotationPairs = ncIntronsPairs, ranks = ncAnnotationRank) %>% 
    mutate(ORFannotation = "intronic")
  
  codingIntronsPairs <- generateTranscriptPairs(transcriptFeature = proteinCodingIntrons, class = "intronic_codingTranscript", annotationTibble = annotationTibble2)
  annotation_codingIntrons <-  prioritizeCodingIntrons(annotationPairs = codingIntronsPairs) %>% 
    mutate(ORFannotation = "intronic")
  
  nonCodingPairs <- generateTranscriptPairs(transcriptFeature = nonCodingExons, class = "non_coding", annotationTibble = annotationTibble2)
  annotation_nonCoding <- prioritizeNonCoding(annotationPairs = nonCodingPairs, ranks = ncAnnotationRank) %>%
    #Annotate as ncRNA, unless it is pseudogenic or antisense
    mutate(ORFannotation = "ncRNA") %>% 
    mutate(ORFannotation = ifelse(transcriptBiotype == "antisense", "antisense", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(transcriptBiotype == "retained_intron", "retained_intron", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(transcriptBiotype == "bidirectional_promoter_lncRNA", "bidirectional_promoter_lncRNA", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(transcriptBiotype == "nonsense_mediated_decay", "nmd", ORFannotation)) %>% 
    mutate(ORFannotation = ifelse(str_detect(transcriptBiotype, 'pseudogene'), "pseudogene", ORFannotation)) 
  
  proteinCodingPairs <- generateTranscriptPairs(transcriptFeature = proteinCodingExons, class = "protein_coding", annotationTibble = annotationTibble2)
  proteinCodingPairsDetailed <- classifyProteinCoding(annotationPairs = proteinCodingPairs, txdb = txdbInput)
  annotation_coding <- prioritizeCoding(annotationPairs = proteinCodingPairsDetailed, ranks = proteinCodingAnnotationRank)
  
  
  annotationCombined <- annotationTibble2 %>% 
    filter(transcriptClass == "intergenic") %>% 
    bind_rows(annotation_ncIntrons) %>% 
    bind_rows(annotation_codingIntrons) %>% 
    bind_rows(annotation_nonCoding) %>% 
    bind_rows(annotation_coding)
  
  return(annotationCombined)
}

#Get IDs of all nORFs in frame with a coding exon
generateInFrameIDs <- function(gffFile, bed6File, txdb) {
  gff0 = read_tsv(gffFile, comment = "#", col_names = F, col_types = 'ccccccccc') 
  gff <- gff0 %>% 
    filter(X3 == "exon") %>% 
    separate(X9, sep = ";", c("a", "b","c","d","e","f","g","h")) %>% 
    dplyr::select(b,e,d, X7) 
  
  colnames(gff) <- c("exon_ID", "start_phase", "end_phase", "strand")
  gff <- gff %>% 
    mutate(exon_ID = str_remove(exon_ID, "Name=")) %>% 
    mutate(end_phase = str_remove(end_phase, "ensembl_end_phase=")) %>% 
    mutate(start_phase = str_remove(start_phase, "ensembl_phase=")) 
  
  bed0 = read_tsv(bed6File, col_names = F, col_types = 'ciiccc') 
  bedFirstExonPlus <-  bed0 %>% 
    filter(X6 == "+") %>% 
    arrange(X1, X2) %>% 
    distinct(X1, X4, .keep_all = T)
  bedFirstExonMinus <- bed0 %>% 
    filter(X6 == "-") %>% 
    arrange(X1, -X3) %>% 
    distinct(X1,X4, .keep_all = T)
  bedFirstExon <- bedFirstExonPlus %>% 
    bind_rows(bedFirstExonMinus) %>% 
    arrange(X1, X2)
  write_tsv(bedFirstExon, "all_38_firstExon.6.bed", col_names = F)
  
  exonsAll <- exons(txdb, use.names=T)
  
  novelORFsExons <- import.bed("all_38_firstExon.6.bed")
  novelORFsExons <- renameSeqlevels(novelORFsExons, gsub("chr","", seqlevels(novelORFsExons)))
  
  inFrame <- findOverlaps(novelORFsExons, exonsAll, type = "any")
  
  inFrameTibble <- tibble(norf_ID = novelORFsExons@elementMetadata@listData$name[inFrame@from],
                          norfStart = novelORFsExons@ranges@start[inFrame@from],
                          norfEnd = novelORFsExons@ranges@start[inFrame@from] + novelORFsExons@ranges@width[inFrame@from] - 1, 
                          exon_ID = exonsAll@ranges@NAMES[inFrame@to],
                          exonStart = exonsAll@ranges@start[inFrame@to],
                          exonEnd =  exonsAll@ranges@start[inFrame@to] + exonsAll@ranges@width[inFrame@to] - 1)  %>% 
    left_join(gff, by = "exon_ID") %>% 
    filter(end_phase != "-1" | start_phase != "-1") 
  
  #Add reading frames for + strand
  inFrameNorfsPlus <- inFrameTibble %>% 
    filter(strand == "+") %>% 
    mutate(norfFrame = norfStart %% 3) %>%
    #If exon has start_phase -1 then shift based on end phase
    #Else shift based on start phase
    mutate(exonFrame = ifelse(start_phase == "-1", (exonEnd - as.integer(end_phase) + 1) %% 3, (exonStart - as.integer(start_phase)) %% 3)) %>% 
    filter(exonFrame == norfFrame) %>% 
    distinct(norf_ID)
  
  inFrameNorfsMinus <- inFrameTibble %>% 
    filter(strand == "-") %>% 
    mutate(norfFrame = norfEnd %% 3) %>%
    #If exon has end_phase -1 then shift based on start phase
    #Else shift based on end phase
    mutate(exonFrame = ifelse(start_phase == "-1", (exonStart + as.integer(end_phase) - 1) %% 3, (exonEnd + as.integer(start_phase)) %% 3)) %>% 
    filter(exonFrame == norfFrame) %>% 
    distinct(norf_ID)
  
  inFrameNorf_IDs <- inFrameNorfsPlus %>% 
    bind_rows(inFrameNorfsMinus)
  
  cds <- cds(human)
  cdsOverlap <- findOverlaps(novelORFs, cds)
  cdsNovelORFs <- novelORFs@elementMetadata@listData$V4[unique(cdsOverlap@from)]

  #Choose only inFrameNorfs that actually overlap cds
  inFrameNorf_IDs <- inFrameNorf_IDs %>% 
    filter(norf_ID %in% cdsNovelORFs)
  
  return(inFrameNorf_IDs)
}
