# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Biostrings and pr2database packages
BiocManager::install("Biostrings")
install.packages("devtools")
devtools::install_github("pr2database/pr2database")

# Load the required libraries
library(pr2database)
library(Biostrings)
library(dplyr)

# Fetch the PR2 database
pr2 <- pr2_database()

version <- packageVersion("pr2database")

# Create a DNAStringSet from the sequences in the PR2 database
ref_dna_seqs <- DNAStringSet(pr2$sequence)

# Create taxonomic names by concatenating columns from the PR2 database
taxon_names <- paste0(pr2$pr2_accession, " ", 
                      paste(pr2$domain, pr2$supergroup, pr2$division, 
                            pr2$subdivision, pr2$class, pr2$order, 
                            pr2$family, pr2$genus, pr2$species, 
                            sep = ";"),";")

# Assign the taxonomic names as the names of the DNAStringSet
names(ref_dna_seqs) <- taxon_names

# Construct filename with version
fasta_filename <- paste0("pr2_", version, "_sequences.fasta")

# Write the sequences to a FASTA file
writeXStringSet(ref_dna_seqs, fasta_filename)

# Filter the PR2 table for rows where 'subdivision' contains 'Foraminifera'
foraminifera_table <- pr2 %>%
  filter(grepl("Foraminifera", subdivision))

# Convert the filtered table to a data frame
foraminifera_table <- as.data.frame(foraminifera_table)

# View the first few rows of the filtered data (optional)
head(foraminifera_table)
