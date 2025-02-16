# Load the DADA2 processing pipeline function from an external script
source("dada2_run.R")  # Ensures modularity and reusability

# Set working directory path and input
path <- getwd()  # Gets the current working directory
t2s_file <- file.path(path, "t2s.csv") # Metadata file mapping tags to samples


# Define file paths for output data
# Create output directory
today_date <- format(Sys.time(), "%Y%m%d_%H%M")
output_dir <- paste0("DADA2_", today_date)
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

asv_table <- file.path(output_dir, "ASV_table.csv") # Output file for ASV table
stats_file <- file.path(output_dir, "filter_stats.csv") # Output file for filtering statistics
repseq_file <- file.path(output_dir, "representative_sequences.fasta") # Output file for representative sequences

# Set processing options
cpus <- TRUE  # Enable multi-threading (set actual number if needed)
pool <- "pseudo"  # DADA2 pooling method for denoising

# Run the DADA2 pipeline with defined inputs and parameters
run_dada2_pipeline(path, t2s_file, asv_table, stats_file, cpus = cpus)
