#accepts the tranformed OTU table (samples x OTUs instead of OTUs x samples as preferred by QIIME)
#returns an OTU table with the prop_cutoff threshold applied
#e.g., a cutoff of 500 would mean any OTU contributing > 500 sequences would be removed
# NOTE: this might be aggressive because of the influence of deleting OTUs on the margin of thresholds
#   i.e., just above the threshold in one treatment, just below the threshold in another treatment
#so, there will often be little change until lower range of this contribution (possibly find a better implementation)
abundance_filter <- function(t_otu_table, otu_abundance_cutoff){
  #calculate the size of an OTU that should be remove (i.e., everything over X)
  #then delete rows that sum to 0
  if(otu_abundance_cutoff <= 1) {return(t_otu_table)} #return if cutoff is going to remove everything

  new_table           <- t_otu_table
  new_table[new_table > otu_abundance_cutoff] <- 0
  excluded_otus       <- rowSums(new_table) == 0
  new_table           <- new_table[rowSums(new_table) != 0,]
  
  #Print updating information to terminal FIX
  #TODO: fix for dropping treatments (as a catch)
  cat("Exclude > ", otu_abundance_cutoff, "\n")
  cat("Dropped treatments: ", paste(row.names(new_table)[excluded_otus]),"\n")
  return(new_table)
}

abundance_filter_tax <- function(pre_otu_table, abundance_cutoff){
  if(abundance_cutoff == 100) {return(pre_otu_table)} #return if cutoff is going to remove nothing
  counts_stop         <- length(pre_otu_table) - 6 #excludes K|P|C|O|F|G columns
  new_table           <- pre_otu_table
  new_table[1:counts_stop][new_table[1:counts_stop] > abundance_cutoff] <- 0 #need to replace over only some columns (maybe convert to matrix then back to data frame?)
  excluded_otus       <- rowSums(new_table[1:counts_stop]) == 0
  #excluded_otus       <- rowSums(subset(new_table, select=1:counts_stop)) == 0
  new_table           <- new_table[rowSums(new_table[1:counts_stop]) != 0,]
  #new_table           <- new_table[rowSums(subset(new_table, select=1:counts_stop)) != 0,]
  
  #Print updating information to terminal FIX
  #TODO: fix for dropping treatments (as a catch)
  cat("Exclude > ", abundance_cutoff, "\n")
  cat("Dropped treatments: ", paste(row.names(new_table)[excluded_otus]),"\n") #seems to be dropping samples not treatments
  return(new_table)
}

###########
#FUNCTIONS#
###########
# Filter sequences that match a model (excluding those that fail to match a SSU model)
ssu_filter <- function(novel_sequences){
  print("Filtering novel sequences against SSU model")
  writeXStringSet(novel_sequences, file = PREFILTER_FASTA)
  filter_command <- paste(SSU_ALIGN, "--no-align", PREFILTER_FASTA, FILTER_DIR, sep = " ")
  system(filter_command, wait = TRUE)
  filtered_sequences <- readRNAStringSet(FILTERED_FASTA)
  return(filtered_sequences)
}

# Align novel sequences to the SSU model (cmalign)
# Merge with reference alignment (esl-alimerge)
align_and_merge <- function(novel_sequences){

  # write novel sequences to novel.fasta
  print("Writing novel sequences to file")
  writeXStringSet(novel_sequences, file = NOVEL_FASTA)
  
  # align novel sequences
  alignment_command <- paste(CMALIGN, "--cpu 2 --mxsize 2056 -o", NOVEL_ALIGNMENT, CMALIGN_MODEL, NOVEL_FASTA, sep = " ")
  print("Calling alignment command")
  print(alignment_command)
  system(alignment_command, wait = TRUE)
  print("Done calling alignment_command")

  # esl-alimerge -o MERGED_ALIGNMENT NOVEL_ALIGNMENT REF_ALIGNMENT
  merge_command <- paste(ESL_ALIMERGE, "--rfonly --outformat afa -o", MERGED_ALIGNMENT, NOVEL_ALIGNMENT, REF_ALIGNMENT, sep = " ")
  system(merge_command)
  #novel_ref_alignment <- readRNAMultipleAlignment(MERGED_ALIGNMENT, format = "stockholm")
  #novel_ref_alignment <- readRNAMultipleAlignment(MERGED_ALIGNMENT)
  novel_ref_alignment <- readRNAStringSet(MERGED_ALIGNMENT)
  return(novel_ref_alignment)
}

make_tree <- function(novel_ref_alignment){
  #dna_alignment <- rna2dna(novel_ref_alignment)
  #writeXStringSet(novel_ref_alignment, seqtype = "DNA", file = TMP_FASTA)
  tree_command <- paste(FASTTREE, "-nt -gtr", MERGED_ALIGNMENT, ">", REF_NOVEL_TREE, sep = " ")
  system(tree_command)
  novel_ref_tree <- read.tree(REF_NOVEL_TREE)
  return(novel_ref_tree)
}

add_to_reference_tree <- function(reference_tree, alignment){
  print("Implement <add_to_reference_tree>")
  # add to reference tree using pplacer
  # use resulting placement and alignment to infer branch lengths using FastTree
  return(1)
}

#### PHYLOGENETIC FUNCTIONS ####

#loosly adapted from: http://ib.berkeley.edu/courses/ib200b/scripts/_R_tree_functions_v1.R
get_parent <- function(tree, node_num)
{
  matching_edges = tree$edge[,2] == node_num
  parent_node_num = tree$edge[,1][matching_edges][1]
  return(parent_node_num)
}

#returns the mrca node of a set of taxa (implement and possibly modify find_novel_clades to return only
# the novel clades - and I could then use this function to return the clade with closest reference)
get_parent_group <- function(tree, list_of_taxa)
{
  #starting_taxon <- get_tip_number(tree, list_of_taxa[1])
  #while still going{
  #  get parent of starting point
  #  if all sequences in group are accounted for{return clade defined by parent of the current node}
  #}
  print("Need to implement <get_parent_group()>")
  return(1)
}
#find the number of a terminal node by taxon name
get_tip_number <- function(tree, tip_label)
{
  return(match(tip_label, tree$tip.label))
}

# returns a list of novel clades (phylo objects) based on the list of novel OTUs
# i.e., what is the phylogenetic distribution of novel OTUs
# NOTE: some novel OTUs are hit multiple times if they are part of the non-novel sister clade of a novel OTU
# TODO: check that all novel OTUs are covered by the novel clades (some reference sequences making it through)
find_novel_clades <- function(tree, novel_otus){
  #checked_novels <- as.character(length(novel_otus))
  checked_novels <- as.character()
  novel_otus_in_clade <- list()
  novel_clade_list <- NULL
  novel_clade_w_refs_list <- NULL
  clade_has_one_terminal <- TRUE
  for(i in 1:length(novel_otus)){
    #print(base::setdiff(checked_novels, novel_otus))
    current_clade <- NULL
    clade_still_novel <- TRUE
    novel_label <- novel_otus[i]
    if(novel_label %in% checked_novels){
      next
    }
    while(clade_still_novel){
      if(is.null(current_clade)){
        clade_has_one_terminal <- TRUE
        novel_node_num <- get_tip_number(tree, novel_label)
        parent_of_novel_clade <- get_parent(tree, novel_node_num)
        current_clade <- extract.clade(tree, parent_of_novel_clade)
        previous_clade <- current_clade
        terminal_nodes <- current_clade$tip.label
        #terminals_ref_clade <- terminal_nodes[!terminal_nodes == novel_label]
        #ref_to_keep <- terminals_ref_clade[1] # keep one leaf in clade sister to novel OTU (to keep parent node of novel for distance calculations)
        #drop_labels <- terminal_nodes[!terminal_nodes %in% c(novel_label, ]
        #previous_clade <- drop.tip(previous_clade, )
      } else{
        clade_has_one_terminal <- FALSE
        previous_clade <- current_clade
        parent_of_novel_clade <- get_parent(tree, parent_of_novel_clade) #node number now representing the crown of the novel clade
        current_clade <- extract.clade(tree, parent_of_novel_clade)
        terminal_nodes <- current_clade$tip.label
      }
      if(length(setdiff(terminal_nodes,novel_otus)) > 0){ # if terminal nodes contain reference OTUs
        clade_still_novel <- FALSE
        novel_clade_list[[length(novel_clade_list) + 1]] <- previous_clade
        novel_clade_w_refs_list[[length(novel_clade_w_refs_list) + 1]] <- current_clade
        #checked_novels <- c(checked_novels, terminal_nodes) #here's the error
        #fixing?:
        if(clade_has_one_terminal){
          checked_novels <- c(checked_novels, novel_label)
          novel_otus_in_clade[[length(novel_otus_in_clade) + 1]] <- list(novel_label)
        }else{
          previous_clade_labels <- previous_clade$tip.label
          #print(previous_clade_labels)
          checked_novels <- c(checked_novels, previous_clade_labels)
          novel_otus_in_clade[[length(novel_otus_in_clade) + 1]] <- list(previous_clade_labels)
        }
      }
    }
  }
  #print(paste(length(novel_clade_list), length(novel_otus_in_clade)), sep = " ")
  new_list <- list("novels"=novel_clade_list, "novels_w_ref"=novel_clade_w_refs_list, "novel_otus_in_clade" = novel_otus_in_clade)
  return(new_list)
}

#calculate the average distance from novel nodes to mrca with ref sequences
#returns NA if tree with references is too large (i.e., novel clade is sister to a huge portion of the tree)
#   - this is a result of dist.nodes failing on too large a tree (overflow)
novelty_distance <- function(clade_with_refs, list_of_novel_otus){
  #list_of_novel_otus <- clade_without_refs$tip.label
  print(paste("Number of terminal nodes in clade:", length(list_of_novel_otus), sep = " "))
  #find the mrca of the novel and reference taxa (node in clade with no edge entering)
  novel_ref_parent <- setdiff(clade_with_refs$edge[,1], clade_with_refs$edge[,2])
  novel_tip_numbers <- lapply(list_of_novel_otus, get_tip_number, tree = clade_with_refs)
  #otu_distances <- dist.nodes(clade_with_refs)[novel_ref_parent, unlist(novel_tip_numbers)] # unlist b/c lapply returns a list, but a vector required
  otu_distances <- numeric(length(list_of_novel_otus))
  otu_distances <- try(dist.nodes(clade_with_refs)[novel_ref_parent, unlist(novel_tip_numbers)], silent=TRUE)
  return(mean(otu_distances))
}

# returns a pssm or scoring model for the novel clade
# to help identify potentially realted contigs in metagenomes
get_pssm <- function(list_of_clade_otus, alignment){
  print("Need to implement <get_pssm()>")
  return(1)
  }

# search metagenome contigs for novel SSU tag
# report visualization of contigs (much more involved analysis here)
tag_metagenome_search <- function(ssu_pssm){
  print("Need to implement the metagenome contig search")
  return(1)
}

#sorts clades for novelty as a function of phylogenetic distance and abundance
sort_clades <- function(clades_list){
  print("Need to implement <sort_clades()>")
  return(1)
}
#### MISC FUNCTIONS ####
# returns vector corresponding to sums of abundances corresponding to OTUs in the novel clades
abundances_by_metadata <- function(list_of_otus, otu_table, metadata_mapping, metadata_category = "all"){
  all_sample_lists <- NULL
  list_of_metadata_values <- "all" #default value for plotting different metadata categories in different columns

  if(metadata_category == "all"){
    all_sample_lists <- metadata_mapping[,1]
    abundance_vector <- numeric(length = 1)
  }else{
    #print(metadata_category)
    abundance_vector <- numeric(length = length(list_of_metadata_values))
    list_of_metadata_values <- unique(metadata_mapping[,metadata_category])
    #print(list_of_metadata_values)
    for(i in 1:length(list_of_metadata_values)){
      temp_list <- NULL
      temp_list <- metadata_mapping[metadata_mapping[,metadata_category] == list_of_metadata_values[i], 1]
      #temp_list <- cleaned_metadata_mapping[cleaned_metadata_mapping[,metadata_category] == list_of_metadata_values[i], 1]
      all_sample_lists[[i]] <- temp_list
    }
  }
  
  if(BIOM){
    subset_otu_table <- subset(otu_table, rownames(otu_table) %in% list_of_otus)    
  }else{
    subset_otu_table <- subset(otu_table, otu_table$X.OTU.ID %in% list_of_otus)
  }
  num_rows_in_subset <- length(subset_otu_table[,1])
  print(paste("num_rows_in_subset:", num_rows_in_subset))
  
  #do I need to convert to a matrix? No, it works as is - the matrix was an awkward work-around before casting labels to a character matrix
  #subset_otu_table <- as.matrix(sapply(subset_otu_table, as.numeric))
  #subset_otu_table <- matrix(as.numeric(unlist(subset_otu_table)),nrow=nrow(subset_otu_table))
  if(metadata_category == "all"){
    # old version: abundance_vector <- sum(subset_otu_table[,-c(1, length(subset_otu_table))])
	abundance_vector <- sum(subset_otu_table)
	print(abundance_vector)
  }else{
    for(i in 1:length(list_of_metadata_values)){
      #sum columns corresponding to the sample list (all_sample_lists[[i]])
      if(num_rows_in_subset == 1){
        temp <- sum(subset_otu_table[,as.character(all_sample_lists[[i]])])
      } else{
        column_sums <- colSums(subset_otu_table[,as.character(all_sample_lists[[i]])])
        temp <- sum(column_sums)
      }
      #abundance_vector[i] <- sum(temp)
      abundance_vector[i] <- temp
    }
  }
  return_list <- list("abundances" = abundance_vector, "metadata_values" = list_of_metadata_values)
  return(return_list)
}

# returns vector corresponding to sums of abundances corresponding to OTUs in the novel clades - BIOM version
abundances_by_metadata_BIOM <- function(clade, otu_table, metadata_mapping, metadata_category = "all"){
  all_sample_lists <- NULL
  
  if(metadata_category == "all"){
    all_sample_lists <- metadata_mapping[,1]
    abundance_vector <- numeric(length = 1)
  }else{
    print(metadata_category)
    list_of_metadata_values <- unique(metadata_mapping[,metadata_category])
    print(list_of_metadata_values)
    for(i in 1:length(list_of_metadata_values)){
      temp_list <- NULL
      temp_list <- metadata_mapping[metadata_mapping[,metadata_category] == list_of_metadata_values[i], 1]
      #temp_list <- cleaned_metadata_mapping[cleaned_metadata_mapping[,metadata_category] == list_of_metadata_values[i], 1]
      all_sample_lists[[i]] <- temp_list
    }
  }
  
  otus_in_clade <- clade$tip.label
  print(otus_in_clade)
  
  subset_otu_table <- subset(otu_table, otu_table$X.OTU.ID %in% otus_in_clade)
  num_rows_in_subset <- length(subset_otu_table[,1])
  print(num_rows_in_subset)
  
  #do I need to convert to a matrix? No, it works as is - the matrix was an awkward work-around before casting labels to a character matrix
  #subset_otu_table <- as.matrix(sapply(subset_otu_table, as.numeric))
  #subset_otu_table <- matrix(as.numeric(unlist(subset_otu_table)),nrow=nrow(subset_otu_table))
  
  if(metadata_category == "all"){
    abundance_vector <- sum(subset_otu_table[,-c(1, length(subset_otu_table))])
  }else{
    for(i in 1:length(list_of_metadata_values)){
      #sum columns corresponding to the sample list (all_sample_lists[[i]])
      if(num_rows_in_subset == 1){
        temp <- sum(subset_otu_table[,as.character(all_sample_lists[[i]])])
      } else{
        column_sums <- colSums(subset_otu_table[,as.character(all_sample_lists[[i]])])
        temp <- sum(column_sums)
      }
      #abundance_vector[i] <- sum(temp)
      abundance_vector[i] <- temp
    }
  }
  return(abundance_vector)
}

# plot novel clade with heatmap of abundances across treatments for each novel node
