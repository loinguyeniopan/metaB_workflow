# metaB_DADA2

# demutiplex
demultRun() demultiplex a list of illumina amplicon libraries (dual-indexed) in parallel.
> head(run_params)
# A tibble: 6 × 4
  library_name primers                          R1                          R2                         
  <chr>        <chr>                            <chr>                       <chr>                      
1 L1_1417      foram37F41_primer14-17_tag.fasta L1_1417_S30_R1_001.fastq.gz L1_1417_S30_R2_001.fastq.gz
2 L2_1417      foram37F41_primer14-17_tag.fasta L2_1417_S31_R1_001.fastq.gz L2_1417_S31_R2_001.fastq.gz

# dada2_workflow
This script is a wrapper for running a DADA2 pipeline on environmental DNA (eDNA) sequencing data in libraries (with "pseodo" pool). 
It is designed to process 
- sequencing reads
- filter and trim
- denoise the sequences
- merge paired-end reads
- remove chimeras
- generate an amplicon sequence variant (ASV) table.

# Key Steps in the Script
# 1. Load the DADA2 pipeline function
The script sources (source("dada2_run.R")) an external R script (dada2_run.R) that contains the actual DADA2 processing pipeline.
This modular approach ensures that the main script remains clean and easy to maintain.

# 2. Define Paths for Input and Output Files
t2sENI.csv: Metadata file that maps sequencing runs to sample names.
>   head(t2s)
  run sample forward reverse
1 L01  A0101    F1-A    17-K
2 L01  A0201    F1-A    17-Q
3 L01  A0301    F1-B    17-P
4 L02  A0102    F1-A    17-K
5 L02  A0202    F1-A    17-Q
6 L02  A0302    F1-B    17-P

ENI_ASV_table.csv: Output table containing ASV counts per sample.
ENI_filter_stats.csv: Output statistics file summarizing read filtering results.

# 3. Set Processing Parameters
cpus: Enables multi-threading for faster processing.
pool: Uses pseudo-pooling during denoising to improve ASV detection.

# 4. Run the DADA2 Pipeline
The script calls run_dada2_pipeline(), passing the required parameters.
This function (from run_dada2.R) performs all processing steps, including:
Quality filtering and trimming
Error model learning
ASV inference
Merging of paired-end reads
Chimera removal
Generation of ASV tables.

