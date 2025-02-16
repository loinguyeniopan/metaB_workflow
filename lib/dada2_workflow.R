# Load the DADA2 processing pipeline function from an external script
source("dada2_run.R")  # Ensures modularity and reusability

# Set working directory path
path <- getwd()  # Gets the current working directory

# Define file paths for input and output data
t2s_file <- file.path(path, "t2sENI.csv")         # Metadata file mapping tags to samples
asv_table <- file.path(path, "ENI_ASV_table.csv") # Output file for ASV table
stats_file <- file.path(path, "ENI_filter_stats.csv") # Output file for filtering statistics

# Set processing options
cpus <- TRUE  # Enable multi-threading (set actual number if needed)
pool <- "pseudo"  # DADA2 pooling method for denoising
lib_list <- unique(t2s$run)  # Identify unique sequencing libraries in metadata

# Display message about the number of libraries being processed
message(paste("DADA2 analysis for", length(lib_list), "libraries."))

# Run the DADA2 pipeline with defined inputs and parameters
run_dada2_pipeline(path, t2s_file, asv_table, stats_file, cpus = cpus)
