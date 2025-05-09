library(dada2)
library(seqinr)
library(Biostrings)
library(ShortRead)
library(R.utils)

# Set the directory containing .gz files
fastq_dir <- getwd()   # adjust if your folder name is different
files <- list.files(fastq_dir, pattern = "\\.gz$", full.names = TRUE)

# Unzip each file
for (f in files) {
  gunzip(f, remove = TRUE, overwrite = TRUE)  # set remove = TRUE to delete .gz after unzip
}


# Set working directory path and input
path <- getwd()  # Gets the current working directory
list.files(path)

# primer check 
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
samples <- sub("_\\d+\\.fastq$", "", basename(fnFs))


names(fnFs) <- samples
names(fnRs) <- samples
run <- "PRJNA1040471"
fwd_primer <- "AAGGGCACCACAAGAACGC"
forward <- "s14F1"

rev_primer <- "GGTCACGTTCGTTGC"
reverse <- "s17"



fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))


# This function counts number of reads in which the primer is found
count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

count_primers(fwd_primer, fnFs[[1]])

count_primers(rev_primer, fnRs[[1]])

cutadapt <- "/Users/ngoc-loinguyen/miniconda3/envs/cutadapt-env/bin/cutadapt"
system2(cutadapt, args = "--version")


# Create an output directory to store the clipped files
cut_dir <- file.path(path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(fnFs))
rev_cut <- file.path(cut_dir, basename(fnRs))

names(fnFs) <- samples
names(fnRs) <- samples

# It's good practice to keep some log files so let's create some
# file names that we can use for those 
cut_logs <- path.expand(file.path(cut_dir, paste0(samples, ".log")))

cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
                   "-G", rev_primer, "-A", fwd_primer_rev,
                   "-n", 2, "--discard-untrimmed")

# Loop over the list of files, running cutadapt on each file.  If you don't have a vector of sample names or 
# don't want to keep the log files you can set stdout = "" to output to the console or stdout = NULL to discard
for (i in seq_along(fnFs)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   fnFs[i], fnRs[i]),
          stdout = cut_logs[i])  
}

# quick check that we got something
head(list.files(cut_dir))

##################################

# Set working directory path and input
path <- cut_dir  # Gets the after cut working directory
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample <- sub("_\\d+\\.fastq", "", basename(fnFs))

# make the t2s file for the run
t2s <- as.data.frame(sample)
t2s$run <- run
t2s <- t2s[, c("run", setdiff(names(t2s), "run"))]
t2s$forward <- forward
t2s$reverse <- reverse
t2s$fwd_primer <- fwd_primer
t2s$rev_primer <- rev_primer

head(t2s)

write.csv(t2s, file = paste0(run, "_t2s.csv"), row.names = FALSE)

t2s_file <- list.files(pattern = "_t2s\\.csv$", full.names = TRUE)


# Read tag-to-sample file
t2s <- read.table(t2s_file, header=TRUE, sep=",")
t2s_name <- sub(".csv", "", basename(t2s_file))


# Define paths
path.filt <- file.path(path, "filterAndTrimed")
dir.create(path.filt, showWarnings = FALSE)


# Define file names after filtering
filtFs <- file.path(path.filt, basename(fnFs))
filtRs <- file.path(path.filt, basename(fnRs))

# Filter and trim
filtering <- filterAndTrim(fwd=fnFs, rev=fnRs,
                           filt=filtFs, filt.rev=filtRs,
                           multithread=TRUE, verbose=TRUE)
filter_stats <- filtering
head(filter_stats)

# Remove samples with zero reads
noReads <- sub("_1.fastq$", "", rownames(filtering)[filtering[, "reads.out"] == 0])
withReads  <- sub("_1.fastq$", "", rownames(filtering)[filtering[, "reads.out"] > 0])
# Tracking no_reads
if (nrow(filtering) != length(t2s$sample))
{
  noReads_track <- FALSE
} else  {noReads_track <- TRUE}

# Keep only valid samples
filtFs_kept <- filtFs[grepl(paste0("^", withReads, collapse = "|"), basename(filtFs))]
filtRs_kept <- filtRs[grepl(paste0("^", withReads, collapse = "|"), basename(filtRs))]
t2s_keep <- t2s[!paste0(t2s_name, "_", t2s$run, "_", t2s$sample) %in% noReads,]

# Learn error models and denoise
lib_sample <- FALSE
by_lib <- TRUE
lib_list <- unique(t2s$run)
pool <- "pseudo"

cpt <- 1

for (i in lib_list)
{
  message(paste("DADA2: learning ", paste(cpt, "/", length(lib_list))))
  # get the same rows from t2s
  if (by_lib)  assign(paste("t2s_keep_",i, sep=""), t2s_keep[grep(paste0("^", i, "$"), t2s_keep$run),])
  
  t2s_tmp <- t2s_keep[grep(paste0("^", i, "$"), t2s_keep$run),]
  assign(paste0("t2s_keep_", i), t2s_tmp)
  
  # Extract sample names
  sample_names <- t2s_tmp$sample
  
  # Fetch files matching any sample name
  filtFs_matched <- filtFs_kept[grepl(paste(sample_names, collapse = "|"), filtFs_kept)]
  filtRs_matched <- filtRs_kept[grepl(paste(sample_names, collapse = "|"), filtRs_kept)]
  
  assign(paste0("filtFs_", i), filtFs_matched)
  assign(paste0("filtRs_", i), filtRs_matched)
  
  
  # to ensure reproducibility
  set.seed(100)
  assign(paste("errFWD_",i, sep=""), learnErrors(get(paste0("filtFs_",i)), nbases = 1e8, multithread=TRUE, randomize=TRUE))
  set.seed(100)
  assign(paste("errREV_",i, sep=""), learnErrors(get(paste0("filtRs_",i)), nbases = 1e8, multithread=TRUE, randomize=TRUE))
  message("DADA2: model ready...")
  message("DADA2: denoising step")
  
  # get the samples
  filtFs_n <- get(paste("filtFs_",i, sep=""))
  filtRs_n <- get(paste("filtRs_",i, sep=""))
  # and the errors profiles
  errF <- get(paste("errFWD_",i, sep=""))
  errR <- get(paste("errREV_",i, sep=""))
  ## derep and dada
  drpF <- derepFastq(filtFs_n)
  drpR <- derepFastq(filtRs_n)
  # rename with sample names
  
  for (j in 1:length(filtFs_n))
  {
    name_full <- filtFs_n[j]
    prefix <- paste(path, "/filterAndTrimed/", sub(".csv", "", tail(strsplit(t2s_file, "/")[[1]],1)), "_",
                    as.character(get(paste0("t2s_keep_",i))[j,"run"]), "_", sep="")
    suffix <- "_fwd_noPrimers.fastq"
    suffix2 <- "_fwd.fastq"
    name <- sub(prefix, "", name_full)
    name <- sub(suffix, "", name)
    name <- sub(suffix2, "", name)
    if (by_lib) names(drpF)[j] <- name
    if (by_lib) names(drpR)[j] <- name
  }
  
  message("DADA2: denoising forward reads")
  ddF <- dada(drpF, err=errF, selfConsist=F, multithread=TRUE, pool = pool)
  message("DADA2: denoising reverse reads")
  ddR <- dada(drpR, err=errR, selfConsist=F, multithread=TRUE, pool = pool)
  merger <- mergePairs(ddF, drpF, ddR, drpR, trimOverhang=TRUE)
  # export merger file
  saveRDS(merger, file.path(path, paste0(i, "_merger_lib.rds")))
  
  message(paste(j, "/", length(filtFs_n), " sample ", name, " is done for : ", i, sep=""))
  cpt <- cpt + 1
}



### Create sequence table, and remove chimeras
all_rds <- sort(list.files(path, pattern="_merger_lib.rds", full.names = TRUE))
all_rds_file <- all_rds_names <- all_rds

for(i in 1:length(all_rds))
{
  all_rds_file[i] <- gsub(paste(path, "/", sep=""), "", all_rds_names[i])
  all_rds_names[i] <- gsub("_merger_lib.rds", "", all_rds_file[i])
}

# If by library (not sample), reconstruct mergers individually and combine
if (by_lib & !lib_sample) {
  all_merger <- c()
  for (i in seq_along(all_rds)) {
    file_to_read <- file.path(path, all_rds_file[i])
    merger_object <- readRDS(file_to_read)
    var_name <- paste0("merger_", all_rds_names[i])
    assign(var_name, merger_object)
    all_merger <- c(all_merger, get(var_name))
  }
} else {
  # Otherwise, read into a named list
  all_merger <- vector("list", length(all_rds_file))
  names(all_merger) <- all_rds_names
  for (i in all_rds_names) {
    file_to_read <- file.path(path, paste0(i, "_merger_lib.rds"))
    all_merger[[i]] <- readRDS(file_to_read)
  }
}

# Make table, chimeras not yet removed
ASV_table <- makeSequenceTable(all_merger)


# count total reads per sample (with chimeras)
withChim <- rowSums(ASV_table)


# Remove the path from the names
names(withChim) <- gsub(paste0("^", path, "/filterAndTrimed/|_1\\.fastq$"), "", names(withChim))


# Consensus chimera removal, recommended
ASV_table_consensus <- removeBimeraDenovo(ASV_table, method = "consensus", multithread=TRUE)

# Now extract the ASV sequences and paste them to make a fasta
ASV_headers <- paste0("ASV", c(1:ncol(ASV_table_consensus)))
seqs <- getSequences(ASV_table_consensus)

library(seqinr)
write.fasta(as.list(seqs), ASV_headers, "repseq_file.fasta")

# now paste the ASV instead of sequences in the matrix before writing it, sort it to match the t2s and transpose
colnames(ASV_table_consensus) <- ASV_headers

rownames(ASV_table_consensus) <- sub("_\\d+\\.fastq$", "", basename(rownames(ASV_table_consensus)))



# create empty vector for samples with no reads
if (length(noReads) > 0 && noReads_track == TRUE)
{
  tmp <- c()
  for (i in 1:length(noReads))
  {
    noReads[i] <- sub(paste0(t2s_name, "_"), "", noReads[i])
    noReads[i] <- sub("_fwd.fastq", "", noReads[i])
    # remove the lib ID from noReads
    for (j in unique(t2s_keep$run)) noReads[i] <- sub(paste0(j, "_"), "", noReads[i])
    tmp <- rbind(tmp, rep(0, length(ASV_headers)))
  }
  rownames(tmp) <- noReads
  # and for the filtering stats
  emptySamples <- c(rep(0,length(noReads)))
  names(emptySamples) <- noReads
  withChim <- c(withChim, emptySamples)
}


# add the samples with no reads in the ASV table
if (length(noReads) > 0  && noReads_track == TRUE) ASV_table_consensus <- rbind(ASV_table_consensus, tmp)

# transpose table and sort table and noChimera count as in the t2s
if (noReads_track) samples_names <- as.character(t2s$sample)
if (!noReads_track) samples_names <- as.character(t2s_keep$sample[paste0(t2s_name, "_", t2s_keep$run, "_", t2s_keep$sample) %in% withReads])

ASV_table_consensus <- t(ASV_table_consensus[samples_names,])

head(rownames(ASV_table_consensus))

ASV_table_consensus <- data.frame(cbind(ASV_ID = ASV_headers, ASV_table_consensus))

# finally write the ASV sorted table
write.table(ASV_table_consensus, file = "asv_table.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")

withChim <- withChim[samples_names]

# adding the statistics on dada2 discarding reads
# cleaning the rownames of filter_stat
for (i in 1:length(rownames(filter_stats)))
{
  rownames(filter_stats)[i] <- sub(paste0(t2s_name, "_"), "", rownames(filter_stats)[i])
  rownames(filter_stats)[i] <- sub("_fwd_noPrimers.fastq", "", rownames(filter_stats)[i])
  rownames(filter_stats)[i] <- sub("_fwd.fastq", "", rownames(filter_stats)[i])
}

for (i in unique(t2s$run)) rownames(filter_stats) <- sub(paste0(i,"_"), "", rownames(filter_stats))

rownames(filter_stats) <- sub("_1\\.fastq$", "", rownames(filter_stats))
# sorting as in the t2s
filter_stats <- filter_stats[samples_names,]
filter_stats <- cbind(sample_ID = rownames(filter_stats), filter_stats, reads.dada2 = withChim)

tmp <- data.frame(ASV_table_consensus[,c(2:ncol(ASV_table_consensus))])
# annoying no numeric
for (i in 1:ncol(tmp)) tmp[,i] <- as.numeric(as.character(tmp[,i]))
filter_stats <- cbind(filter_stats, reads.dada2.noChimera = colSums(tmp))

write.table(filter_stats, file = "stats.tsv", quote = F, sep="\t", row.names = F, fileEncoding = "UTF-8")
