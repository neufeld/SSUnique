#!/bin/bash

echo "Installing base dependencies for SSUnique...\n"

# SOFTWARE:
# FastTree - http://www.microbesonline.org/fasttree/FastTree
# infernal - http://selab.janelia.org/software/infernal/infernal-1.1.1.tar.gz 
# ssualign - ftp://selab.janelia.org/pub/software/ssu-align/ssu-align-0.1.tar.gz
# base R - sudo apt-get install r-base r-base-dev
# HMMER: - http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz

# LIBRARIES:
# phyloseq

# DATA:
# Living Tree Project: http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_123/LTPs123_SSU.compressed.fasta

mkdir -p software
mkdir -p downloads
mkdir -p data

cd downloads
# Installing FastTree
echo "Installing FastTree"
if hash FastTree 2>/dev/null; then
    echo "FastTree already installed"
else
    echo "Installing from http://www.microbesonline.org/fasttree/FastTree"
    wget "http://www.microbesonline.org/fasttree/FastTree"
    cp FastTree /usr/local/bin/
    cd ..
fi

# Installing Infernal
echo "Installing Infernal"
if hash cmalign 2>/dev/null; then
    echo "Infernal already installed"
else
    echo "Installing from http://selab.janelia.org/software/infernal/infernal-1.1.1.tar.gz"
    wget "http://selab.janelia.org/software/infernal/infernal-1.1.1.tar.gz"
    tar -xv infernal-1.1.1.tar.gz
    cd infernal-1.1.1
    ./configure
    make
    make install
    cd ..
fi

# Installing ssu-align
echo "Installing ssu-align"
if hash ssu-align 2>/dev/null; then
    echo "ssu-align already installed"
else
    echo "Installing from ftp://selab.janelia.org/pub/software/ssu-align/ssu-align-0.1.tar.gz"
    wget "ftp://selab.janelia.org/pub/software/ssu-align/ssu-align-0.1.tar.gz"
    tar -xv ssu-align0.1.tar.gz
    cd ssu-align-0.1
    ./configure
    make
    make install
    cd ..
fi

# Installing HMMER
echo "Installing HMMER"
if hash hmmer 2>/dev/null; then
    echo "HMMER already installled"
else
    echo "Installing from http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz"
    wget "http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz"
    tar -cvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
    cd hmmer-3.1b2-linux-intel-x86_64
    ./configure
    make
    make install
    cd ..
fi

# Installing base R

# Installing OligoArrayAux (for DECIPHER)
echo "Instalilng OligoArrayAux (DECIPHER dependency)"
if hash hybrid-min 2>/dev/null; then
    echo "OligoArrayAux installed"
else
    echo "Installing from http://unafold.rna.albany.edu/cgi-bin/OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz"
    wget "http://unafold.rna.albany.edu/cgi-bin/OligoArrayAux-download.cgi?oligoarrayaux-3.8.tar.gz"
    tar -cvf oligoarrayaux-3.8.tar.gz
    cd oligoarrayaux-3.8
    ./configure
    make
    make install
    cd ..
fi

# Getting and prepping standards data
# download Living Tree Project data
echo "Setting up LTP SSU reference data"
cd data
wget "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_123/LTPs123_SSU.compressed.fasta"

# remove dashes and realign using Infernal

echo -e "The base dependencies of SSUnique have been installed and is ready to run. A detailed manual and sample data can be found at: \"github.com/mdjlynch/SSUnique\"
