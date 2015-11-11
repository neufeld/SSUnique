# Introduction
SSUnique is a pipeline for the identification of phylogenetic novelty in 16S microbiome studies. It has been successfully used to identify significant phylogentic novelty in Earth Microbiome Project and Human Microbiome Project data (Lynch and Neufeld in review). SSUnique has been tested in Linux and Mac environments.

A reference manual and sample workflow are also available.

# Dependencies
SSUnique relies on several pieces of software and reference data for basic functionality:

- FastTree
- Infernal
- ssu-align
- HMMER
- R (including BioStrings, phyloseq, and ggplot2)
- Reference data (Living Tree Project)

# Installation
SSUnique dependencies and reference data can be installed/checked manually or using the install script:

./base_check.sh

# References
Lynch et al. 2012
Lynch and Neufeld in review
