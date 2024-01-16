#!/usr/bin/env Rscript

# Run interactively for now

if (!require('yaml')) install.packages('yaml'); library('yaml')

#args <- commandArgs(trailingOnly=TRUE)
#if(length(args) == 0){
#        cat("Please provide a YAML file. \n")
#        cat("Usage: Rscript --vanilla scr00-00-replicate-check.R sample_library.yaml \n")
#        stop()
#}

#params <- read_yaml(file.path(args[1]))
params <- read_yaml(file.path("meta", "sample_library.yaml")) # For testing keep it commented in production

samples <- list.dirs(file.path(params$demult_dir), recursive = FALSE, full.names = FALSE)

split_list <- strsplit(samples, "_")

df_samples <- do.call(rbind, lapply(split_list, function(x) data.frame("subject_id" = x[1], 
                                                                       "visit" = x[2], 
                                                                       "replicate" = x[3])))
grouped_counts <- aggregate(x = df_samples$subject_id, by = list(df_samples$subject_id), FUN = length)
colnames(grouped_counts) <- c("subject_id", "replicate_count")

df_out <- merge(df_samples, grouped_counts)

write.table(df_out, file.path("results", "sample-replicate-check.tsv"), sep = "\t",
row.names = FALSE, quote= FALSE)
