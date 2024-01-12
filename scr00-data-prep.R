#!/usr/bin/env Rscript

fastq_dir <- file.path("data/AACFV5MHV/EVD68_PhIP_Seq/")
sample_id <- read.csv2(file.path("data/PhIPseq_EVD68_Sample_Sheet_121323.txt"), 
                       sep = "\t", header = TRUE, check.names = FALSE)
sample_id <- sample_id[["Sample_ID"]]
sample_id <- sample_id[sample_id != ""]


data_dir <- file.path("data", "demux-folders")
dir.create(data_dir, recursive=TRUE, showWarnings=FALSE)

for(i in 1:length(sample_id)){
        samp <- sample_id[i]
        cat("Copying files for:",i, ":" , samp, "\n")
        fastq_files <- sort(list.files(fastq_dir, pattern = paste0(samp, "*"), 
                                       full.names = TRUE, recursive = FALSE))
        samp1 <- file.path(data_dir, paste0(samp, "_1"))
        samp2 <- file.path(data_dir, paste0(samp, "_2"))

        dir.create(samp1, recursive = TRUE, showWarnings=FALSE)
        dir.create( samp2, recursive = TRUE, showWarnings=FALSE)

        file.copy(fastq_files[1], samp1, overwrite = TRUE)
        file.copy(fastq_files[2], samp1, overwrite = TRUE)
        
        file.copy(fastq_files[3], samp2, overwrite = TRUE)
        file.copy(fastq_files[4], samp2, overwrite = TRUE)
}






