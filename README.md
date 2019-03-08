# nORF-data-prep
The following is a command-line walkthrough for processing [OpenProt](http://www.openprot.org/) ([Brunet et al., 2018](https://doi.org/10.1093/nar/gky936)) and [sorfs.org](http://www.sorfs.org/) ([Olexiouk et al., 2018](https://doi.org/10.1093/nar/gkx1130)) human novel open reading frame (nORF) databases into lab usable formats.

## Prerequisites
1. Clone this repository
```
git clone https://github.com/PrabakaranGroup/nORF-data-prep.git
```
2. R with the following packages installed:
* tidyverse
* Biostrings

## Processing OpenProt

#### Step 1: Download Files

The code below downloads files with the following parameters from OpenProt (http://www.openprot.org/p/download):
* Release: 1.3
* Species: Homo Sapiens
* Assembly: GRCg38.p5
* Protein Type: AltProts, Isoforms and Refprots
* Annotation: Ensembl (GRCh38.83) + RefSeq (GRCh38.p7) + UniProt (2017-09-27)

Download, unzip, and rename the 'All predicted' `.tsv`, `.bed`, `.dna.fasta`, and `.fasta` files:
```
#TSV file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.tsv.zip
unzip human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.tsv.zip
mv human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.tsv openProtAllPredicted.tsv

#BED file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.bed.zip
unzip human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.bed.zip
mv human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.bed openProtAllPredicted_38.bed

#DNA FASTA file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.dna.fasta.zip
unzip human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.dna.fasta.zip
mv human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7.dna.fasta openProtAllPredicted_38.dna.fasta

#Protein FASTA file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.fasta.zip
unzip human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.fasta.zip
mv human-openprot-r1_3-refprots+altprots+isoforms-grch38.83+grch38.p7+uniprot2017_09_27.fasta openProtAllPredicted_38.fasta

#Clean up
rm human-openprot-r1_3-refprots+altprots+isoforms*
```
#### Step 2: Process Files

Use the provided `processOpenProt.R` to create custom bed and FASTA files. 
Specifically, this script will:
 * A) Extract OpenProt entries labelled as 'AltProt' with mass spec AND/OR ribo-seq evidence. These are the entries that qualify as nORFs.
 * B) Output subsetted `.gtf` and `.bed` files
 * C) Optionally, if passed `2peptide` argument, will create the same files subsetted to the entries having 2+ MS peptides of evidence

```
#OpenProt with any mass-spec/ribo-seq evidence
Rscript processOpenProt.R

#OpenProt with 2+ peptides of mass-spec evidence
Rscript processOpenProt.R 2peptide
```


## Processing sorfs.org

#### Step 1: Download Files

The code below downloads files with the following parameters from sorfs.org (http://www.sorfs.org/BioMart):
* Database: Homo Sapiens
* Floss Classification: Good, Extreme
* Main Attributes: Sorf ID, Chromosome, Sorf Start, Sorf End, Strand, Spliced Start Parts, Spliced Stop Parts, Start Codon, Sorf Length, AA sequence, Transcript sequence, Biotype, Annotation, Ensembl Transcript ID
```
wget -O sorfsDownload.txt 'http://biomart.biobix.be/martservice/results?query=<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="BioMart" config="Human"><Filter name="human__classification_104" value="Extreme,Good" filter_list=""/><Attribute name="human__sorf_id_104"/><Attribute name="human__chr_104"/><Attribute name="human__sorf_begin_104"/><Attribute name="human__sorf_end_104"/><Attribute name="human__strand_104"/><Attribute name="human__start_parts_104"/><Attribute name="human__stop_parts_104"/><Attribute name="human__sorf_length_104"/><Attribute name="human__start_codon_104"/><Attribute name="human__aa_seq_104"/><Attribute name="human__biotype_104"/><Attribute name="human__annotation_104"/><Attribute name="human__id_104"/><Attribute name="human__tr_seq_104"/></Dataset></Query>'
```
Optionally, the code below downloads a similar file for entries with peptide evidence. It uses the same attributes as above plus:
* Re_Spin_attributes: PRIDE file
```
wget -O sorfsDownload_MS.txt 'http://biomart.biobix.be/martservice/results?query=<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="BioMart" config="Human"><Filter name="human__classification_104" value="Extreme,Good" filter_list=""/><Attribute name="human__sorf_id_104"/><Attribute name="human__chr_104"/><Attribute name="human__sorf_begin_104"/><Attribute name="human__sorf_end_104"/><Attribute name="human__strand_104"/><Attribute name="human__start_parts_104"/><Attribute name="human__stop_parts_104"/><Attribute name="human__sorf_length_104"/><Attribute name="human__start_codon_104"/><Attribute name="human__aa_seq_104"/><Attribute name="human__tr_seq_104"/><Attribute name="human__biotype_104"/><Attribute name="human__annotation_104"/><Attribute name="human__id_104"/><Attribute name="ReSpin__file_106"/></Dataset></Query>'
```

#### Step 2: Process Files

Use the provided `processSorfs.R` to create custom bed and FASTA files. 
Specifically, this script will:
 * A) Convert genomic coordinates and splicing info into bed12 format
 * B) Remove duplicate entries based on matching genomics coordinates
 * C) Output subsetted `.gtf` and `.bed` files
 * D) Optionally, if passed `massSpec` argument, will create the same files subsetted to the entries with MS evidence
 
```
#sorfs.org all unique entries
Rscript processSorfs.R

#sorfs.org subsetted to entries with mass-spec evidence
Rscript processSorfs.R massSpec
```

## Authors

Written by Matt Neville (Cambridge Department of Genetics), Narendra Meena (Department of Biology, IISER Pune), Chaitanya Erady (Cambridge Department of Genetics), and Robin Kohze (Cambridge Department of Genetics).

## Acknowledgments
Supervisor: Sudhakaran Prabakaran (https://prabakaran-group.org/)
