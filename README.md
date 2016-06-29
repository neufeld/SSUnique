# Introduction
SSUnique is a pipeline for the identification of phylogenetic novelty in 16S microbiome studies. It has been successfully used to identify significant phylogentic novelty in Earth Microbiome Project and Human Microbiome Project data (Lynch and Neufeld in review). SSUnique has been tested in Linux and Mac environments.

A reference manual and sample workflow are also available.

# Dependencies
SSUnique relies on several pieces of software and reference data for basic functionality:

- [FastTree] (http://www.microbesonline.org/fasttree) (tested with v. 2.1.3 SSE3)
- [infernal] (http://eddylab.org/infernal/)  (tested with v. 1.1.1)
- [ssualign] (http://eddylab.org/software/ssu-align/) (tested with v. 0.1.1)
- [base R] (https://www.r-project.org/) (R >= v. 3.2)
- [HMMER] (http://hmmer.org/) (v. 3.1)

- Reference Data [SILVA Living Tree Project] (http://www.arb-silva.de/projects/living-tree/) (tested with v. 119)

# Installation
SSUnique dependencies and reference data can be installed/checked manually or using the install scripts:

## Dependencies

sudo ./base_check.sh

## Reference Data
sudo ./reference_data.sh

# Manual

A brief manual, including a sample workflow, can be found:

# References
1. Lynch et al. 2012. Targeted recovery of novel phylogenetic diversity from next-generation sequence data. ISME J. 6:2067–2077 doi:10.1038/ismej.2012.50
2. Lynch, MDJ and JD Neufeld. Automated exploration of microbiology’s most wanted in global high-throughput sequence data (in review)
