run_dada2_pipeline <- function(path, t2s_file, asv_table, stats_file, cpus = TRUE) {
 
  library(dada2)
  library(seqinr)
  
  message("DADA2: Start filter and trim...")
  
  # Read tag-to-sample file
  t2s <- read.table(t2s_file, header=TRUE, sep=",")
  t2s_name <- sub(".csv", "", basename(t2s_file))
  
  # Define paths
  path.filt <- file.path(path, "filterAndTrimed")
  dir.create(path.filt, showWarnings = FALSE)
  
  # Get FASTQ files
  fwd <- sort(list.files(path, pattern = "_fwd.fastq", full.names=TRUE))
  rev <- sort(list.files(path, pattern = "_rev.fastq", full.names=TRUE))
  
  # Define file names after filtering
  filtFs <- file.path(path.filt, basename(fwd))
  filtRs <- file.path(path.filt, basename(rev))
  
  # Filter and trim
  filtering <- filterAndTrim(fwd=fwd, rev=rev, filt=filtFs, filt.rev=filtRs, multithread=cpus, verbose=TRUE)
  
  # Remove samples with zero reads
  noReads <- sub("_fwd.fastq$", "", rownames(filtering)[filtering[, "reads.out"] == 0])
  withReads  <- sub("_fwd.fastq$", "", rownames(filtering)[filtering[, "reads.out"] > 0])
  # Tracking no_reads
  if (nrow(filtering) != length(t2s$sample))
  {
    noReads_track <- FALSE
  } else  {noReads_track <- TRUE}
  
  # Keep only valid samples
  filtFs_kept <- filtFs[!grepl(paste0("^", noReads, collapse = "|"), basename(filtFs))]
  filtRs_kept <- filtRs[!grepl(paste0("^", noReads, collapse = "|"), basename(filtRs))]
  t2s_keep <- t2s[!paste0(t2s_name, "_", t2s$run, "_", t2s$sample) %in% noReads,]
  
  message("DADA2: filterAndTrim done.")
  
  # Learn error models and denoise
  lib_list <- unique(t2s$run)
  cpt <- 1
  for (i in lib_list) {
    message(paste("DADA2: learning library", cpt,"/", length(lib_list)))
    
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
    
    message("\n DADA2: denoising forward reads")
    ddF <- dada(drpF, err=errF, selfConsist=F, multithread=cpus, pool = pool)
    message("\n DADA2: denoising reverse reads")
    ddR <- dada(drpR, err=errR, selfConsist=F, multithread=cpus, pool = pool)
    merger <- mergePairs(ddF, drpF, ddR, drpR, trimOverhang=TRUE)
    
    # export merger file
    
    saveRDS(merger, file.path(path, paste0(i, "_merger_lib.rds")))
    for (k in 1:length(merger)) saveRDS(merger[k], file.path(path, paste0(names(merger)[k], "_merger.rds")))
    
    message(paste("\n Total ", j, "/", length(filtFs_n), " samples is done for : ", i, sep=""))
    
    cpt <- cpt + 1
    
  }
  
  
  ### Create sequence table, and remove chimeras
  all_rds <- sort(list.files(path, pattern="_merger.rds", full.names = TRUE))
  all_rds_file <- all_rds_names <- all_rds
  
  # reconstruct the merger
  for(i in 1:length(all_rds))
  {
    all_rds_file[i] <- gsub(paste(path, "/", sep=""), "", all_rds_names[i])
    all_rds_names[i] <- gsub("_merger.rds", "", all_rds_file[i])
  }
  
  all_merger <- c()
  for(i in 1:length(all_rds))
  {
    assign(paste0("merger_",all_rds_names[i]), readRDS(file.path( paste0(all_rds_file[i]))))
    all_merger <- c(all_merger, get(paste0("merger_",all_rds_names[i])))
  }
  
  # Make table, chimeras not yet removed
  ASV_table <- makeSequenceTable(all_merger)
  
  # count total reads per sample (with chimeras)
  withChim <- rowSums(ASV_table)
  
  # Remove chimeras
  ASV_table_consensus <- removeBimeraDenovo(ASV_table, method="consensus", multithread=TRUE)
  
  # Now extract the ASV sequences and paste them to make a fasta
  ASV_headers <- paste0("ASV", c(1:ncol(ASV_table_consensus)))
  seqs <- getSequences(ASV_table_consensus)
  write.fasta(as.list(seqs), ASV_headers, repseq_file)
  
  
  # now paste the ASV instead of sequences in the matrix before writing it, sort it to match the t2s and transpose
  colnames(ASV_table_consensus) <- ASV_headers
  
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
  withChim <- withChim[samples_names]
  ASV_table_consensus <- data.frame(cbind(ASV_ID = ASV_headers, ASV_table_consensus))
  
  
  # Save ASV table as CSV
  write.csv(ASV_table_consensus, file=asv_table, row.names = FALSE, fileEncoding="UTF-8")
  
  
  # adding the statistics on dada2 discarding reads
  filter_stats <- filtering
  
  # cleaning the rownames of filter_stat
  for (i in 1:length(rownames(filter_stats)))
  {
    rownames(filter_stats)[i] <- sub(paste0(t2s_name, "_"), "", rownames(filter_stats)[i])
    rownames(filter_stats)[i] <- sub("_fwd_noPrimers.fastq", "", rownames(filter_stats)[i])
    rownames(filter_stats)[i] <- sub("_fwd.fastq", "", rownames(filter_stats)[i])
  }
  
  for (i in unique(t2s$run)) rownames(filter_stats) <- sub(paste0(i,"_"), "", rownames(filter_stats))
  # sorting as in the t2s
  filter_stats <- filter_stats[samples_names,]
  filter_stats <- cbind(sample_ID = rownames(filter_stats), filter_stats, reads.dada2 = withChim)
  
  tmp <- data.frame(ASV_table_consensus[,c(2:ncol(ASV_table_consensus))])
  # annoying no numeric
  for (i in 1:ncol(tmp)) tmp[,i] <- as.numeric(as.character(tmp[,i]))
  filter_stats <- cbind(filter_stats, reads.dada2.noChimera = colSums(tmp))
  
  # Save filter_stats as CSV
  write.csv(filter_stats, file = stats_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  message("DADA2 pipeline completed for ", length(lib_list) , " library/libraries.")
}
