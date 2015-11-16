# Setting universal constants, loading required libraries, and sourcing SSUnique scripts

#####################
# REQUIRED PACKAGES #
#####################
library(optparse)
library(ape)
library(biom)
library(phyloseq)
library(ggplot2)
library(Biostrings)
library(reshape) # Needed?

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
METADATA_CAT        <- "all"#"pH_FineRange"#"pH_LargeRange"#"Site_ID"#"Land_Use_Type"#"Depth_cm" #"HMPbodysubsite" #"all" #"ENV_FEATURE"
BIOM                <- "TRUE"

# program paths
FASTTREE            <- '/usr/local/bin/FastTreeMP'
PPLACER             <- 'pplacer'
CMALIGN             <- '/Winnebago/mike/_bioinformatics_tools/cmalign'
ESL_ALIMERGE        <- '/Winnebago/mike/_bioinformatics_tools/esl-alimerge'
RAXML               <- '/home/mdjlynch/bioinformatics_tools/misc_tools/bin/raxmlHPC-PTHREADS-SSE3'
SSU_ALIGN           <- '/usr/local/bin/ssu-align'

# program arguments
RAXML_PARAMS        <- '/home/mdjlynch/Research/PostDoc/working/SSUnique/standards/ssunique_standards/RAxML_binaryModelParameters.LTP90ucPARAMS'

# reference standards
REF_ALIGNMENT       <- "/Winnebago/mike/_working/ssunique/ref_data/LTPs119_SSU.cmalign.stk"
REF_NOVEL_TREE      <- paste(RESULTS_DIR, "reference_plus_novel.tre", sep = "")
PPLACER_REFPKG      <- '/home/mdjlynch/Research/PostDoc/working/SSUnique/standards/ssunique_standards/LTP90ucSSUnique.refpkg'
CMALIGN_MODEL       <- '/Winnebago/mike/_working/ssunique/ref_data/bacteria-0p1.infernal.cm'
EPA_TREE            <- 'RAxML_labelledTree.ref_plus_novel'
CLEAN_EPA_TREE      <- 'RAxML_labelledTree.ref_plus_novel.newick' #likely don't need

#SSUnique outputs
PREFILTER_FASTA     <- paste(RESULTS_DIR, "prefiltered.fasta", sep = "")
FILTERED_FASTA      <- paste(FILTER_DIR, "filter_results.bacteria.fa", sep = "")
TMP_FASTA           <- paste(RESULTS_DIR, "tmp.fasta", sep = "")
ALIGN_OUT           <- paste(RESULTS_DIR, "tmp.stk", sep = "")
#MERGED_ALIGNMENT    <- paste(RESULTS_DIR, "merged.stk", sep = "")
MERGED_ALIGNMENT    <- paste(RESULTS_DIR, "merged.afa", sep = "")
NOVEL_FASTA         <- paste(RESULTS_DIR, "novel.fasta", sep = "")
NOVEL_ALIGNMENT     <- paste(RESULTS_DIR, "novel.stk", sep = "")
TMP_TREE            <- paste(RESULTS_DIR, "tmp.tre", sep = "")
REF_NOVEL_SUFFIX    <- "ref_plus_novel"
CLADES_OUTPUT       <- paste(RESULTS_DIR, "clades.txt", sep = "")
HIGHLIGHT_TREE      <- paste(RESULTS_DIR, "highlight_tree.xml", sep = "")
