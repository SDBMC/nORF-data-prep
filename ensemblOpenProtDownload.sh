#!/bin/bash

#TSV file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.83.tsv.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.83.tsv.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.83.tsv dataFiles/ensemblOpenProtAllPredicted.tsv

#BED file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.83.bed.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.83.bed.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.83.bed dataFiles/ensemblOpenProtAllPredicted_38.bed

#Protein FASTA file
wget http://www.openprot.org/download/files/1.3/human-openprot-r1_3-altprots+isoforms-grch38.83.fasta.zip
unzip human-openprot-r1_3-altprots+isoforms-grch38.83.fasta.zip
mv human-openprot-r1_3-altprots+isoforms-grch38.83.fasta dataFiles/ensemblOpenProtAllPredicted_38.fasta

#Clean up
rm human-openprot-r1_3-altprots+isoforms-grch38.83.*
