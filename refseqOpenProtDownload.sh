#!/bin/bash

#TSV file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.p7.tsv.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.p7.tsv.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.p7.tsv refseqOpenProtAllPredicted.tsv

#BED file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.p7.bed.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.p7.bed.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.p7.bed refseqOpenProtAllPredicted_38.bed

#Protein FASTA file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.p7.fasta.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.p7.fasta.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.p7.fasta refseqOpenProtAllPredicted_38.fasta

#Clean up
rm human-openprot-r1_3-altprots+isoforms-grch38.p7.*

