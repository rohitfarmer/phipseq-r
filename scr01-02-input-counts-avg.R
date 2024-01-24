#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr01-03-input-counts-avg.R input_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "input_library.yaml")) # For testing keep it commented in production

# Read count files
count_files <- list.files(file.path(params$output_count_dir), pattern = ".tsv",
                          full.names = TRUE)

df_all_counts <- data.frame()
for(i in 1:length(count_files)){
        f <- count_files[i]
        coln <-  tools::file_path_sans_ext(basename(f))
        df_counts <- read.csv2(f, sep = "\t", header = TRUE, check.names = FALSE,
                               col.names = c("id", coln))
        if(i == 1){
                df_all_counts <- df_counts
        }else{
                df_all_counts <- merge(df_all_counts, df_counts)
        }
}

df_all_counts$Input <- as.integer(round(rowMeans(df_all_counts[, 2:3])))

df_out <- df_all_counts[, c("id", "Input")]

out_file <- file.path(params$output_count_dir, "input-counts-avg.tsv")
write.table(df_out, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
