#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')
if (!require('doMC')) install.packages('doMC'); library('doMC')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr02-01-normalize.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "sample_library_mystat_thru.yaml")) # For testing keep it commented in production

# Functions

user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}

# Define input folders
output_merged_count_dir <- params$output_merged_count_dir

# Define output folders
output_normalized_counts_dir <- params$output_normalized_counts_dir

# Create folders

if(dir.exists(output_normalized_counts_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_normalized_counts_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_normalized_counts_dir, "\n")
        dir.create(output_normalized_counts_dir, recursive = TRUE, showWarnings = FALSE)
}


# Read files
files <- list.files(path = output_merged_count_dir, full.names = TRUE, pattern = "*.tsv")
registerDoMC(detectCores(all.tests = FALSE, logical = FALSE) - 1)
null_out <- foreach(file_i = 1:length(files)) %dopar% {
        input <- files[file_i]
        cat("Processing:", basename(input), "\n")
        input_dat <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

        col1 <- input_dat[[2]] / (sum(input_dat[[2]]) / 1e6)
        col1 <- round(col1)

        col2 <- input_dat[[3]] / (sum(input_dat[[3]]) / 1e6)
        col2 <- round(col2)

        output_dat <- data.frame(input_dat[[1]], col1, col2)
        colnames(output_dat) <- colnames(input_dat)

        output_file_name <- file.path(output_normalized_counts_dir, paste0(tools::file_path_sans_ext(basename(input)), "_normalized.tsv"))
        write.table(output_dat, file = output_file_name, 
                    sep = "\t", row.names = FALSE, quote = FALSE)
}



