# TODO: convert system calls to system2 and direct output to log file. Note, cannot append, so use temporary files
# TODO: add maximum or mean length of sequences in clade to the ranking algorithm
# TODO: when outputting the ranked clades, delete novel OTU terminal nodes in the known seed clade

# Main script for SSUnique, controlling program flow
source("ssunique_settings.R")

# parse command line arguments
option_list <- list(
                    make_option(c("-b", "--biom_table"), action = "store", default = "otu_table.biom",
                                help = "biom-formatted OTU table [default: otu_table.biom]"),
                    make_option(c("-m", "--mapping"), action = "store", default = "mapping.txt",
                                help = "QIIME-formatted metadata mapping file [default: mapping.txt]"),
                    #make_option(c("-c", "--category"), action = "store", default = "all",
                    #            help = "mapping category for analysis [default: all]"),
                    make_option(c("-s", "--sequences"), action = "store", default = "rep_set.fna",
                                help = "sequences for OTU representative set [default = rep_set.fna]"),
                    make_option(c("-x", "--exclude_emp"), action = "store", default = "env:",
                                help = "string for subsetting EMP data (not generally required)")
                    )

opts <- parse_args(OptionParser(option_list = option_list))

OTU_TABLE <- opts$biom_table
METADATA  <- opts$mapping
SEQUENCES <- opts$sequences
# TODO: add category here instead of through ssunique_settings.R print(opts$category)
EMP_EXCLUDE <- opts$exclude_emp # TODO: not generally needed

# create results directories
dir.create("ssunique_results", showWarnings = FALSE)
dir.create(file.path("ssunique_results", "ranked"), showWarnings = FALSE)

# Start log file by coppying ssunique_settings.R
file.copy(from = "ssunique_settings.R", to = LOG_FILE)

# send SSUnique call to log file
args <- commandArgs(trailingOnly = TRUE)
#write(paste("SSUnique v.", VERSION, sep = ""), file = LOG_FILE)
#write(args, file = LOG_FILE, append = TRUE)

# Start log file by coppying ssunique_settings.R
file.copy(from = "ssunique_settings.R", to = LOG_FILE)

# Write information to log file:
#write(paste("otu table:", OTU_TABLE, "\nmetadata:", METADATA, "\nsequences:", SEQUENCES, 
#            sep = " "), file = LOG_FILE, append = TRUE)
#write(paste("Global variables:\nabundance threshold:", ABUNDANCE_THRESHOLD, "\nreference alignment:", 
#            REF_ALIGNMENT, "\nalignment model:", CMALIGN_MODEL, sep = " "), file = LOG_FILE, append = TRUE)

# read OTU table and mapping file
print("Loading/processing biom file and metadata")
otu_data <- import_biom(OTU_TABLE)
print("done loading biom file")
metadata_mapping  <- read.table(METADATA, header = TRUE, sep = "\t", comment.char = "", quote = "\"")
print("done loading mapping file")
# TODO: is this subsetting required? Not generally, but for EMP subsets, yes.
#metadata_mapping <- metadata_mapping[which(tolower(metadata_mapping$ENV_BIOME) == "envo:tropical and subtropical moist broadleaf forest biome"),]
metadata_mapping <- metadata_mapping[which(tolower(metadata_mapping$ENV_BIOME) == EMP_EXCLUDE),]

# table of total abundances for each unique metadata feature studied
metadata_abundances <- sum_table_by_metadata(otu_table(otu_data), METADATA_CAT)

# subset OTU table to novel OTUs (novel at order, which is a superset of novelty at higher ranks)
# NOTE: NA is the "unclassified"
#rename columns since EMP10000 doesn't read in KPCOFGS
colnames(tax_table(otu_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#colnames(tax_table(otu_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
novels <- subset_taxa(otu_data, is.na(Class))
abundant_novels <- prune_taxa(taxa_sums(novels) >= ABUNDANCE_THRESHOLD, novels)
abundant_novel_otus <- taxa_names(abundant_novels)

# free up some memory
remove(otu_data)

print("Loading sequence data")
sequences <- readDNAStringSet(SEQUENCES)
novel_sequences <- sequences[abundant_novel_otus]
# use load_sequence_subset, but with an added error check for OTUs that don't have rep_set representatives
# free up some memory
remove(sequences)
# filter sequences for non-SSU artifacts (in response to novel clades hitting phiX and other scaffolds)
# using ssu-align for filtering
filtered_sequences <- ssu_filter(novel_sequences)
filtered_otus <- names(filtered_sequences)

# Skip time intensive alignment and tree building if files already present
if(file.exists(REF_NOVEL_TREE)){
  novel_plus_ref_aln <- readRNAStringSet(MERGED_ALIGNMENT)
  novel_plus_ref_tree <- read.tree(REF_NOVEL_TREE)
} else{
  # read in sequences, filter by novel OTUs
  #print("Loading sequence data")
  #sequences <- readDNAStringSet(SEQUENCES)
  #novel_sequences <- sequences[abundant_novel_otus]
  # use load_sequence_subset, but with an added error check for OTUs that don't have rep_set representatives
  # free up some memory
  #remove(sequences)
  # filter sequences for non-SSU artifacts (in response to novel clades hitting phiX and other scaffolds)
  # using ssu-align for filtering
  #filtered_sequences <- ssu_filter(novel_sequences)
  #filtered_otus <- names(filtered_sequences)
  
  # align novel sequences to model, merge with reference alignment
  print("Aligning novel sequences to model")
  novel_plus_ref_aln <- align_and_merge(filtered_sequences)
  print("Building tree")
  novel_plus_ref_tree <- make_tree(novel_plus_ref_aln)  
}

# find novel clades
# novel_clades contains novels, novels_w_ref, novel_otus containing novel tree, novel tree with references
#   and novel otus (the novel tree might get deleted since issues arise when the novel clade should only be a leaf node)
print("Finding novel clades")
novel_clades <- find_novel_clades(novel_plus_ref_tree, filtered_otus)
print(paste(length(novel_clades$novels_w_ref), "novel clades observed", sep = " "))

# calculate distances of novel clades
# http://stackoverflow.com/questions/13832221/how-to-be-alerted-about-the-ongoing-progress-of-a-loop-lapply
print("Calculating clade distances")
#distances <- mapply(novelty_distance, clade_with_refs = novel_clades$novels_w_ref, clade_without_refs = novel_clades$novels)
# alternatively, with a for loop
distances <- numeric(length(novel_clades$novels_w_ref))
for(i in 1:length(distances)){
  if(i %% 100 == 0){print(i)}
  distances[i] <- novelty_distance(clade_with_refs = novel_clades$novels_w_ref[[i]], list_of_novel_otus = unlist(novel_clades$novel_otus_in_clade[[i]]))
}

# rank clade novelty (initial ranking done by f(n) = pdist x abundance)
# calculate the ranking factor for each clade, sort that list,
#   giving an association of rank and unranked position for each clade,
#   which gives an easy accessor to a ranking of novelty
ranking <- numeric(length(distances))
abundances_for_ranking <- numeric(length(distances))
ranking_value <- 0

# plot clade novelty
print("Parse clade novelty by metadata category")
abundance_data <- NULL
# temp variable before sorting out plotting - abundance_data table cannot handle metadata subsets
# probably can be removed now
metadata_set <- "all"

# for(i in 1:length(novel_clades$novels_w_ref)){
#   #plot(novel_clades$novels[[i]])
#   abund_by_meta <- abundances_by_metadata(unlist(novel_clades$novel_otus_in_clade[[i]]), otu_table(abundant_novels), metadata_mapping, METADATA_CAT)
#   print(paste("abundances:", abund_by_meta$abundances))
#   temp_abundance_data <- cbind.data.frame(abund_by_meta$metadata_values, abund_by_meta$abundances, 
#                                           rep(distances[i], length(abund_by_meta$metadata_values)))
#   abundance_data <- rbind(abundance_data, temp_abundance_data)
#   
#   abundances_for_ranking[i] <- sum(abund_by_meta$abundances)
# }

# modify so that abundance_data now has a column for proportions
for(i in 1:length(novel_clades$novels_w_ref)){
  #plot(novel_clades$novels[[i]])
  abund_by_meta <- abundances_by_metadata(unlist(novel_clades$novel_otus_in_clade[[i]]), otu_table(abundant_novels), metadata_mapping, METADATA_CAT)
  print(paste("abundances:", abund_by_meta$abundances))
  temp_abundance_data <- cbind.data.frame(abund_by_meta$metadata_values, abund_by_meta$abundances, 
                                          rep(distances[i], length(abund_by_meta$metadata_values)))
  abundance_data <- rbind(abundance_data, temp_abundance_data)
  
  abundances_for_ranking[i] <- sum(abund_by_meta$abundances)
}

colnames(abundance_data) <- c("metadata", "abundance", "pdist")

proportions <- numeric(length(abundance_data[,1]))
for(i in 1:length(abundance_data[,1])) {
  #proportions[i] <- abundance_data$abundance[i]/metadata_abundances$TABLE_SUM[match(metadata_abundances$METADATA, "all")]
  proportions[i] <- abundance_data$abundance[i]/metadata_abundances$TABLE_SUM[metadata_abundances$METADATA == abundance_data$metadata[i]]
}

# add proportions to the abundance_data table (gets output to file)
abundance_data <- cbind(abundance_data, proportions)

# novelty plotting
plot_novelty_sized(abundance_data, METADATA_CAT)
plot_novelty_facet(abundance_data, METADATA_CAT)
# FIX PROPORTIONAL PLOTS
plot_novelty_sized_proportion(abundance_data, METADATA_CAT, metadata_abundances)
plot_novelty_facet_proportion(abundance_data, METADATA_CAT, metadata_abundances)

# novelty ranking by descending order
#   convert each value to a z-score and then sum values, doubling the weight for phylogenetic distance
# possibly try another ranking scheme (potentially just the average of the pdist rank and abundance rank)
###################
# Z-score ranking #
###################
# scaled_abundances <- scale(abundances_for_ranking)
# scaled_distances <- scale(distances)
# ranking_variable <- (scaled_distances * 2) + scaled_abundances
# ranking_table <- cbind(ranking_variable, 1:length(ranking_variable))
# colnames(ranking_table) <- c("RANKING", "ORIGINAL_POSITION")
# ranking_table <- data.frame(ranking_table)
# ranking_table <- ranking_table[with(ranking_table,order(-RANKING)),]

###########################
# Sorted averaged ranking #
###########################
# TODO: check that the scaling of 10 is appropriate here
distance_ranking <- rank(-distances) * 10
#distance_ranking <- distance_ranking * 2
abundance_ranking <- rank(-abundances_for_ranking)
avg_ranking <- rowMeans(cbind(distance_ranking, abundance_ranking))
ranking_variable <- rank(-avg_ranking)
ranking_table <- cbind(ranking_variable, 1:length(distances))
colnames(ranking_table) <- c("RANKING","ORIGINAL_POSITION")
ranking_table <- data.frame(ranking_table)
ranking_table <- ranking_table[with(ranking_table,order(-RANKING)),]

# output clade novelty and abundance data to file
if(METADATA_CAT == "all") {
  filename <- paste(RESULTS_DIR, METADATA_CAT, "_abundance_novelty.table", sep = "")
  # change metadata from "all" to the biom filename
  temp_split <- strsplit(OTU_TABLE, "/")[[1]]
  all_label <- paste(METADATA_CAT, tail(temp_split, n = 1), sep = ":")
  modified_abundance <- abundance_data
  modified_abundance[,1] <- rep(all_label, length(modified_abundance[,1]))
  write.table(modified_abundance, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
  rm(modified_abundance)
} else {
  filename <- paste(RESULTS_DIR, METADATA_CAT, "_abundance_novelty.table", sep = "")
  write.table(abundance_data, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}

# output novel clades to text file
for(i in 1:length(novel_clades$novels_w_ref)){
  filename <- paste(RESULTS_DIR, "trees_novels_w_refs.tre", sep = "")
  write.tree(novel_clades$novels_w_ref[[i]], file = filename)
}

# output treefile (with refs) and alignment subset (without refs) of novel clades
clade_novelty_rank <- 0
num_left_fill <- nchar(length(distances))
for(i in 1:length(ranking_table[,1])){
  clade_novelty_rank <- ranking_table$ORIGINAL_POSITION[i]
  # write treefile (with reference sequences)
  filled_rank <- paste("ranked/ranked_novel_clade", formatC(i, width = num_left_fill, flag = "0"), sep = "_")
  filename <- paste(RESULTS_DIR, filled_rank, "_pdist_", distances[clade_novelty_rank], "_w_ref.tre", sep = "")
  write.tree(novel_clades$novels_w_ref[[clade_novelty_rank]], file = filename)
  
  # write alignment (without reference sequences)
  #otus_in_clade <- novel_clades$novels[[clade_novelty_rank]]$tip.label
  #novel_plus_ref_aln contains aligned novel and reference sequences
  clade_sequences <- novel_plus_ref_aln[unlist(novel_clades$novel_otus_in_clade[[clade_novelty_rank]])]
  filename <- paste(RESULTS_DIR, filled_rank, "_pdist_", distances[clade_novelty_rank], ".afa", sep = "")
  writeXStringSet(clade_sequences, file = filename)
  # added hmmbuild command
  hmm_model <- paste(RESULTS_DIR, filled_rank, "_pdist_", distances[clade_novelty_rank], ".hmm", sep = "")
  hmmbuild_command <-  paste("hmmbuild ", hmm_model, filename)
  system(hmmbuild_command)

  clade_novelty_rank <- clade_novelty_rank + 1
}

# calculate PRT (permenantly rare taxa), output table to file

# calculate TRT (transiently rare taxa), output table to file

# calculate CRT (conditionally rare taxa), output table to file
