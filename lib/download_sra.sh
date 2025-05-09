#!/bin/bash

# Set project accession
PROJECT=PRJEB35877

# Step 1: Get run info from SRA and save to file
esearch -db sra -query ${PROJECT} | efetch -format runinfo | tee runinfo.csv > ${PROJECT}_runinfo.csv

# Step 2: Extract only the Run column (skip header)
tail -n +2 runinfo.csv | cut -d ',' -f1 | grep -E '^SRR|^ERR|^DRR' > srr_list.txt

# Step 3: Download each .sra file using prefetch
cat srr_list.txt | xargs prefetch

# Step 4: Convert .sra files to .fastq.gz using fastq-dump
cat srr_list.txt | xargs -I{} fastq-dump --split-files --gzip {}
