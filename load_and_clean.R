#load libraries (find a more elegant way to check for libraries)
library(ggplot2)
library(gplots)
library(vegan)
library(ape)
library(gridExtra)

#### LOADING DATA ####
load_biom_file <- function(filename){
  print("Need to implement <load_biom_file()>")
  biom_file <- 1
  return(biom_file)
}

# load sequence file, keeping sequences named in the novel_otus list
# deprecated for now, but add check for sequences that are not in the fasta file and throw an error
load_sequence_subset <- function(seq_filename, subset_names = "ALL"){
  sequences <- readDNAStringSet(filepath = seq_filename)
  print(length(subset_names))
  if(toupper(subset_names) == "ALL"){
    seq_subset <- sequences
  } else{
    #seq_subset <- sequences[sapply(subset_names, function(x) grep(x, names(sequences)))]
  }
  return(seq_subset)
}

# give the sum of the entire table or distributed across metadata categories
# takes the OTU table directly instead of the phyloseq object (e.g., otu_table(phyloseq_object))
sum_table_by_metadata <- function(otu_table_matrix, metadata_cat = "all"){
  if(metadata_cat == "all"){
    full_sum <- sum(otu_table_matrix)
    return_table <- data.frame("all", full_sum, stringsAsFactors = FALSE)
    colnames(return_table) <- c("METADATA", "TABLE_SUM")
  } else {
    list_of_metadata_values <-  unique(metadata_mapping[,metadata_cat])
    metadata_values <- character(length(list_of_metadata_values))
    metadata_sums <- numeric(length(list_of_metadata_values))
    for(i in 1:length(list_of_metadata_values)){
      temp_list <- NULL
      temp_list <- metadata_mapping[metadata_mapping[,metadata_cat] == list_of_metadata_values[i], 1]
      temp_sums <- colSums(otu_table_matrix[,as.character(temp_list)])
      metadata_values[i] <- as.character(list_of_metadata_values[i])
      metadata_sums[i] <- sum(temp_sums)
    }
    return_table <- data.frame(metadata_values, metadata_sums, stringsAsFactors = FALSE)
    #TODO (mdjlynch): this data frame has levels that interfere with downstream calculations (may be solved)
    #return_table <- data.frame(metadata_values, metadata_sums)
    colnames(return_table) <- c("METADATA", "TABLE_SUM")
  }
  return(return_table)
}

#############
# OTU TABLE #
#############

# #stacked plot from previous analysis (as guide):
# #otu_plot <- ggplot(otu.table,aes(Encoding,Count,fill=Phylum),alpha=Genus) + geom_bar(position="fill",stat="identity",colour="white") + coord_flip()
# 
# # Abundance Threshold   Taxonomic Cleveland Plot
# # -----|-----           Order1 *************
# #                       Order2 *********
# # Metadata Category     Order3 ****
# # 1:x 2: 3: 4:          Order4 **
# #
# #
# # IMPLEMENT AS A SHINY CHART OR SHINYAPP EMBEDDED IN A DOCUMENT
# # Create an OTU table with taxonomy split up over the last seven columns, K | P | C | O | F | G | S
# # Note, this section and the next can be better optimized (since both are using similar data processing)
# # Note, keep otu table as data frame, converting to matrix only when doing the ordination
# # process OTU table, splitting taxonomy (this is ugly and hackish, but works right now - fix later)
# # from: http://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame
# otu_table_taxsplit <- within(otu_table, Consensus.Lineage<-data.frame(do.call('rbind', strsplit(as.character(Consensus.Lineage), ';', fixed=TRUE))))
# otu_table_taxsplit <- cbind(otu_table_taxsplit, otu_table_taxsplit$Consensus.Lineage$X2, otu_table_taxsplit$Consensus.Lineage$X3,
#                             otu_table_taxsplit$Consensus.Lineage$X4, otu_table_taxsplit$Consensus.Lineage$X5,
#                             otu_table_taxsplit$Consensus.Lineage$X6, otu_table_taxsplit$Consensus.Lineage$X7)
# 
# # rename columns from Consensus.Lineage.X1-X7 to KINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS
# # Exclude SPECIES for now
# # TODO: All of this hacking and splicing might be better dealt with when reading in the file (change
# # ';' to tabs, write to temp file, then read.table). Check: http://www.r-bloggers.com/using-r-reading-tables-that-need-a-little-cleaning/
# # Works as is for now
# otu_table$Consensus.Lineage.X1 <- NULL # remove the 'Root' column
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X2"] <- "KINGDOM"
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X3"] <- "PHYLUM"
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X4"] <- "CLASS"
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X5"] <- "ORDER"
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X6"] <- "FAMILY"
# names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X7"] <- "GENUS"
# #names(otu_table_taxsplit)[names(otu_table_taxsplit)=="otu_table_taxsplit$Consensus.Lineage$X8"] <- "SPECIES"
# 
# # if cell doesn't start with k__ (or whatever), replace with Unclassified
# otu_table_taxsplit$KINGDOM  <- gsub("^((?!k__).)*$", "Unclassified", otu_table_taxsplit$KINGDOM, perl=TRUE)
# otu_table_taxsplit$PHYLUM   <- gsub("^((?!p__).)*$", "Unclassified", otu_table_taxsplit$PHYLUM, perl=TRUE)
# otu_table_taxsplit$CLASS    <- gsub("^((?!c__).)*$", "Unclassified", otu_table_taxsplit$CLASS, perl=TRUE)
# otu_table_taxsplit$ORDER    <- gsub("^((?!o__).)*$", "Unclassified", otu_table_taxsplit$ORDER, perl=TRUE)
# otu_table_taxsplit$FAMILY   <- gsub("^((?!f__).)*$", "Unclassified", otu_table_taxsplit$FAMILY, perl=TRUE)
# otu_table_taxsplit$GENUS    <- gsub("^((?!g__).)*$", "Unclassified", otu_table_taxsplit$GENUS, perl=TRUE)
# #otu_table_taxsplit$SPECIES  <- gsub("^((?!s__).)*$", "Unclassified", otu_table_taxsplit$SPECIES, perl=TRUE)
# 
# #delete Consensus.Lineage
# otu_table_taxsplit$Consensus.Lineage <- NULL
