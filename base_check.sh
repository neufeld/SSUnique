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
    # finish install here
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
    # finish install here
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
    # finish install here
    cd ..
fi

# Installing base R

# Getting and prepping standards data
# download Living Tree Project data
echo "Setting up LTP SSU reference data"
cd data
wget "http://www.arb-silva.de/fileadmin/silva_databases/living_tree/LTP_release_123/LTPs123_SSU.compressed.fasta"

# remove dashes and realign using Infernal


rm -rf downloads

# MetaAnnotate script:
software=`pwd`/software
PATH="${PATH}:${HOME}/.local/bin"

mkdir -p downloads
mkdir -p software

echo "Installing pip.\n"
if [ ! `which pip` ] ; then
  cd downloads
  wget "https://bootstrap.pypa.io/get-pip.py"
  python get-pip.py --user --ignore-installed
  cd ..
fi

echo "Installing EMBOSS transeq.\n"
if [ ! `which transeq` ] ; then
  cd downloads
  wget "ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz"
  tar -xzf EMBOSS-6.5.7.tar.gz
  cd EMBOSS*
  ./configure --without-x
  make
  cp -R emboss/ "$software"/
  ln -s "$software"/emboss/transeq ~/.local/bin/transeq
  cd ..
  cd ..
fi


echo "Installing python packages through pip.\n"
~/.local/bin/pip install --user numpy --ignore-installed
~/.local/bin/pip install --user celery --ignore-installed
~/.local/bin/pip install --user taxtastic --ignore-installed
~/.local/bin/pip install --user lxml --ignore-installed
~/.local/bin/pip install --user python-gflags --ignore-installed
~/.local/bin/pip install --user ete2 --ignore-installed

echo "Installing KronaTools.\n"
if [ ! `which ktImportText` ] ; then
  cd included_software/KronaTools-2.5/
  ./install.pl --prefix "${HOME}/.local/"
  cd ..
  cd ..
fi


echo "Installing HMMER & Easel mini-applications.\n"
if [ ! `which hmmsearch` ] | [ ! `which esl-sfetch` ] ; then
  cd downloads
  wget "ftp://selab.janelia.org/pub/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz"
  tar -xzf hmmer-3.1b1-linux-intel-x86_64.tar.gz
  cd hmmer*
  ./configure
  make
  cp -R . "$software"/hmmer/
  ln -s "$software"/hmmer/binaries/hmmstat ~/.local/bin/hmmstat
  ln -s "$software"/hmmer/binaries/hmmsearch ~/.local/bin/hmmsearch
  ln -s "$software"/hmmer/binaries/hmmalign ~/.local/bin/hmmalign
  ln -s "$software"/hmmer/binaries/esl-reformat ~/.local/bin/esl-reformat
  ln -s "$software"/hmmer/binaries/esl-sfetch ~/.local/bin/esl-sfetch
  cd ..
  cd ..
fi

echo "Installing USEARCH.\n"
if [ ! `which usearch` ] ; then
  ln -s `pwd`"/included_software/usearch" ~/.local/bin/usearch
fi

echo "Installing FastTreeMP.\n"
if [ ! `which FastTreeMP` ] ; then
  cd downloads
  wget "http://www.microbesonline.org/fasttree/FastTreeMP"
  mv FastTreeMP ~/.local/bin/
  chmod a+x ~/.local/bin/FastTreeMP
  cd ..
fi

echo "Installing pplacer and guppy.\n"
if [ ! `which guppy` ] ; then
  cd downloads
  wget "http://matsen.fhcrc.org/pplacer/builds/pplacer-v1.1-Linux.tar.gz"
  tar -xzf pplacer-v1.1-Linux.tar.gz 
  cd pplacer*
  mv pplacer ~/.local/bin/
  mv guppy ~/.local/bin/
  cd ..
  cd ..
fi

echo "Downloading and indexing taxonomy info.\n"
if [ ! -e data/taxonomy.pickle ] ; then
  cd precompute
  wget "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  tar -zxf taxdump.tar.gz
  grep 'scientific name' names.dmp > trimmed.names.dmp
  python make_taxonomy_pickle.py
  cd ..
fi

echo "Downloading and indexing gi number to taxid mappings.\n"
if [ ! -e data/gi_taxid_prot.dmp ] ; then
  cd precompute
  wget "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz"
  gunzip gi_taxid_prot.dmp.gz
  mv gi_taxid_prot.dmp ../data/
  cd ..
fi

rm -rf downloads
rm -f precompute/gc.prt
rm -f precompute/readme.txt
rm -f precompute/taxdump.tar.gz

echo "$HOME/.local/bin/" > path.txt

echo -e "The base dependencies of SSUnique have been installed and is ready to run. A detailed manual and sample data can be found at: \"github.com/mdjlynch/SSUnique\"
