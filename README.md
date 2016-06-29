# Introduction
SSUnique is a pipeline for the identification of phylogenetic novelty in 16S microbiome studies. It has been successfully used to identify significant phylogentic novelty in Earth Microbiome Project and Human Microbiome Project data (Lynch and Neufeld in review). SSUnique has been tested in Linux and Mac environments.

A reference manual and sample workflow are also available.

# Dependencies
SSUnique relies on several pieces of software and reference data for basic functionality:

- FastTree - http://www.microbesonline.org/fasttree/FastTree
- infernal - http://eddylab.org/infernal/infernal-1.1.1.tar.gz 
- ssualign - http://eddylab.org/software/ssu-align/ssu-align-0.1.1.tar.gz
- base R - sudo apt-get install r-base r-base-dev
- HMMER: - http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz

- Reference data (Living Tree Project)

# Installation
SSUnique dependencies and reference data can be installed/checked manually or using the install script:

sudo ./base_check.sh
./reference_data.sh

# References
Lynch et al. 2012
Lynch and Neufeld in review
