# Setting universal constants, loading required libraries, and sourcing SSUnique scripts

#####################
# REQUIRED PACKAGES #
#####################

# CRAN packages, R >= 3.1 (ggplot2)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(optparse, ape, biom, ggplot2, reshape)

# Bioconductor packages, R >= 3.2 (phyloseq)
# Must have Bioconductor installed, see: https://www.bioconductor.org/install/
pacman::p_load(Biostrings, phyloseq)

#########################
# LOAD SSUnique SOURCES #
#########################
source("align_and_merge.R")
source("load_and_clean.R")
source("novelty_functions.R")
source("novelty_plotting_and_output.R")
source("rare_functions.R")
source("rare_plotting_and_output.R")

####################
# SET UP CONSTANTS #
####################

# general
VERSION             <- "0.1"
RESULTS_DIR         <- "ssunique_results/"
FILTER_DIR          <- paste(RESULTS_DIR, "filter_results/", sep = "")
LOG_FILE            <- paste(RESULTS_DIR, "ssunique.log", sep = "")
CLEANUP             <- FALSE
TAX_RANK            <- "Order"
ABUNDANCE_THRESHOLD <- 1 # exclude OTUs that sum below this threshold - useful for very large datasets
METADATA_CAT        <- "SOURCE" #e.g., "HMPbodysubsite", "all", "ENV_FEATURE"
BIOM                <- "TRUE"

# program paths
FASTTREE            <- 'FastTree'
CMALIGN             <- 'cmalign'
ESL_ALIMERGE        <- 'esl-alimerge'
SSU_ALIGN           <- 'ssu-align'

# reference standards
REF_ALIGNMENT       <- "ref_data/LTPs119_SSU.cmalign.stk"
REF_NOVEL_TREE      <- paste(RESULTS_DIR, "reference_plus_novel.tre", sep = "")
CMALIGN_MODEL       <- 'ref_data/SSU_rRNA_bacteria.cm'

#SSUnique outputs
PREFILTER_FASTA     <- paste(RESULTS_DIR, "prefiltered.fasta", sep = "")
FILTERED_FASTA      <- paste(FILTER_DIR, "filter_results.bacteria.fa", sep = "")
TMP_FASTA           <- paste(RESULTS_DIR, "tmp.fasta", sep = "")
ALIGN_OUT           <- paste(RESULTS_DIR, "tmp.stk", sep = "")
MERGED_ALIGNMENT    <- paste(RESULTS_DIR, "merged.afa", sep = "")
NOVEL_FASTA         <- paste(RESULTS_DIR, "novel.fasta", sep = "")
NOVEL_ALIGNMENT     <- paste(RESULTS_DIR, "novel.stk", sep = "")
TMP_TREE            <- paste(RESULTS_DIR, "tmp.tre", sep = "")
REF_NOVEL_SUFFIX    <- "ref_plus_novel"
CLADES_OUTPUT       <- paste(RESULTS_DIR, "clades.txt", sep = "")
HIGHLIGHT_TREE      <- paste(RESULTS_DIR, "highlight_tree.xml", sep = "")
