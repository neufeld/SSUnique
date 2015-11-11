library(ggplot2)
library(reshape)

# novelty table (used in plot_novelty_sized() and plot_novelty_facet()):

# distances abundance_metadata1 abundance_metadata2 ...

# plot clade novelty along a single axis, with different columns for different metadata categories
# metadata_category is just a variable for the filename
# TODO: figure out a better way to scale the size-based plotting
plot_novelty_sized <- function(novelty_table, metadata_category = "all"){
  filename = paste(RESULTS_DIR, metadata_category, "_clade_diversity_sized.pdf", sep = "")
  pdf(file = filename, paper = "a4r")
  p <- ggplot(novelty_table, aes(x = metadata, y = pdist)) + geom_point(aes(size = abundance), color = "blue", alpha = 0.5) + scale_size(range = c(0,25)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  
  plot(p)
  dev.off()
}

# plot clade novelty as a dot plot, facetted by metadata type
# metadata_category is just a variable for the filename
plot_novelty_facet <- function(novelty_table, metadata_category = "all"){
  #ggplot(novelty_melt, aes(x = distances, y = value)) + geom_point() + facet_grid(. ~ variable)
  filename = paste(RESULTS_DIR, metadata_category, "_clade_diversity_facet.pdf", sep = "")
  pdf(file = filename, paper = "a4r")
  p <- ggplot(novelty_table, aes(x = pdist, y = sqrt(abundance))) + geom_point(color = "blue", alpha = 0.5) + facet_wrap(~metadata)
  plot(p)
  dev.off()
}

plot_novelty_sized_proportion <- function(novelty_table, metadata_category = "all", metadata_sums){
  filename = paste(RESULTS_DIR, metadata_category, "_clade_diversity_proportion.pdf", sep = "")
  prop_novelty_table <- novelty_table
  #add new column transforming abundance to proportion - try speedup with an apply function
  #for(i in 1:length(prop_novelty_table$abundance)){
  #  prop_novelty_table$prop[i] <- prop_novelty_table$abundance[i]/metadata_sums$TABLE_SUM[metadata_sums$METADATA == prop_novelty_table$metadata[i]]
  #}
  pdf(file = filename, paper = "a4r")
  #p <- ggplot(prop_novelty_table, aes(x = metadata, y = pdist)) + geom_point(aes(size = prop), color = "blue", alpha = 0.5) + scale_size(trans = "log10", range = c(0,10)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  p <- ggplot(prop_novelty_table, aes(x = metadata, y = pdist)) + geom_point(aes(size = proportions), color = "blue", alpha = 0.5) + scale_size(range = c(0,20)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  plot(p)
  dev.off()
}

plot_novelty_facet_proportion <- function(novelty_table, metadata_category = "all", metadata_sums) {
  filename = paste(RESULTS_DIR, metadata_category, "_clade_diversity_facet_propotion.pdf", sep = "")
  prop_novelty_table <- novelty_table
  #for(i in 1:length(prop_novelty_table$abundance)) {
  #  prop_novelty_table$prop[i] <- prop_novelty_table$abundance[i]/metadata_sums$TABLE_SUM[metadata_sums$METADATA == prop_novelty_table$metadata[i]]
  #}
  pdf(file = filename, paper = "a4r")
  p <- ggplot(novelty_table, aes(x = pdist, y = proportions)) + geom_point(color = "blue", alpha = 0.5) + facet_wrap(~metadata)
  plot(p)
  dev.off()
}
