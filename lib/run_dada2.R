# Load required packages
library(dada2)
library(seqinr)

# Define paths
path <- getwd()
t2s_file <- file.path(path, "t2sENI.csv")
asv_table <- file.path(path, "ENI_ASV_table.tsv")
stats_file <- file.path(path, "ENI_filter_stats.tsv")

# Set CPU cores
cpus <- TRUE 

# Read tag-to-sample file
t2s <- read.table(t2s_file, header=TRUE, sep=",")
t2s_name <- sub(".csv", "", basename(t2s_file))

# Determine if one library per sample
lib_sample <- length(unique(t2s$run)) == length(unique(t2s$sample))
lib_list <- unique(t2s$run)

# Define filtering paths
path.filt <- file.path(path, "filterAndTrimed")
dir.create(path.filt, showWarnings = FALSE)  # Ensure directory exists

# Get FASTQ files
fwd <- sort(list.files(path, pattern = "_fwd.fastq", full.names=TRUE))
rev <- sort(list.files(path, pattern = "_rev.fastq", full.names=TRUE))

# Define file names after filtering
filtFs <- file.path(path.filt, basename(fwd))
filtRs <- file.path(path.filt, basename(rev))

# Filter and trim
filtering <- filterAndTrim(fwd=fwd, rev=rev, filt=filtFs, filt.rev=filtRs, multithread=cpus, verbose=TRUE)

# Remove samples with zero reads
noReads <- rownames(filtering)[filtering[, "reads.out"] == 0]
withReads <- rownames(filtering)[filtering[, "reads.out"] > 0]

# Keep only valid samples
filtFs_kept <- filtFs[!basename(filtFs) %in% noReads]
filtRs_kept <- filtRs[!basename(filtRs) %in% noReads]

t2s_keep <- t2s[!paste0(t2s_name, "_", t2s$run, "_", t2s$sample) %in% gsub("_fwd\\.fastq$", "", noReads), ]

message("DADA2: filterAndTrim done.")




# DADA2 processing
message("DADA2: Learning error model(s) step")

cpt <- 1

for (i in lib_list)
{
  message(paste("DADA2: learning ", paste(cpt, "/", length(lib_list))))
  # fetch here either the list of sample to be processed, or the unique sample from the list
  assign(paste("filtFs_",i, sep=""), filtFs_kept[grep(paste0(i, "_"), filtFs_kept)])
  assign(paste("filtRs_",i, sep=""), filtRs_kept[grep(paste0(i, "_"), filtRs_kept)])
 
   # get the same rows from t2s
  assign(paste("t2s_keep_",i, sep=""), t2s_keep[grep(paste0("^", i, "$"), t2s_keep$run),])

  # to ensure reproducibility
  set.seed(100)
  assign(paste("errFWD_",i, sep=""), learnErrors(get(paste0("filtFs_",i)), nbases = 1e8, multithread=cpus, randomize=TRUE))
  set.seed(100)
  assign(paste("errREV_",i, sep=""), learnErrors(get(paste0("filtRs_",i)), nbases = 1e8, multithread=cpus, randomize=TRUE))
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
   names(drpF)[j] <- name
   names(drpR)[j] <- name
  }
  
  message("DADA2: denoising forward reads")
  ddF <- dada(drpF, err=errF, selfConsist=F, multithread=cpus, pool="pseudo")
  message("DADA2: denoising reverse reads")
  ddR <- dada(drpR, err=errR, selfConsist=F, multithread=cpus, pool="pseudo")
  merger <- mergePairs(ddF, drpF, ddR, drpR, trimOverhang=TRUE)
  
  # export merger file

saveRDS(merger, file.path(path, paste0(i, "_merger.rds")))
  
  cpt <- cpt + 1
}

