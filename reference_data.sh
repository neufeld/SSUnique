#!/bin/bash

echo "Setting up data dependencies for SSUnique...\n"

# DATA:
# Living Tree Project: http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_123/LTPs123_SSU.compressed.fasta

mkdir -p ref_data

# get cm model
cd ref_data
wget "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.1/Rfam.cm.gz"
cmfetch Rfam.cm.gz SSU_rRNA_bacteria > SSU_rRNA_bacteria.cm
cmfetch Rfam.cm.gz SSU_rRNA_archaea > SSU_rRNA_archaea.cm
cmfetch Rfam.cm.gz > SSU_rRNA_eukarya.cm SSU_rRNA_eukarya > SSU_rRNA_eukarya.cm
rm Rfam.cm.gz
cd ..

# Getting and prepping standards data
# download Living Tree Project data
echo "Setting up LTP SSU reference data"
cd ref_data
wget "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_123/LTPs123_SSU.compressed.fasta"

# remove dashes and realign using Infernal
cat LTPs123_SSU.compressed.fasta | tr "\t" "_" | tr -d " " > temp.fasta
cmalign -o LTPs119_SSU.cmalign.stk SSU_rRNA_bacteria.cm temp.fasta
rm LTPs123_SSU.compressed.fasta
rm temp.fasta
cd ..

echo "Finished setting up data dependencies for SSUnique\n"
