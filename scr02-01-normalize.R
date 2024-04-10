#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')
if (!require('doMC')) install.packages('doMC'); library('doMC')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr02-filter-merge.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "ntc_library.yaml")) # For testing keep it commented in production

# Functions

filter_reads <- function(output_count_dir, output_passed_count_dir){

        # Listing files in the Sample Counts directory
        files <- list.files(path = file.path(output_count_dir), full.names = TRUE)

        # Looping through files
        cores <-  detectCores()
        registerDoMC(cores)
        bh <- foreach (i =1:length(files)) %dopar% {
                f <- files[i]
                if (!grepl("^\\.", basename(f))) { # Excluding hidden files (similar to startswith('.') in Python)
                        name <- tools::file_path_sans_ext(basename(f))
                        cat("Sample:",name, "\n")

                        # Reading the tab-separated file
                        sample_counts <- read.csv2(f, sep = "\t", header = TRUE, check.names = FALSE)

                        # Calculating MIN_READS
                        min_reads <- (nrow(sample_counts) * 10) / 2

                        # Summing aligned counts (column-wise)
                        aligned_counts <- colSums(sample_counts[2], na.rm = TRUE)
                        #cat(aligned_counts[1], "\n")

                        # Saving the file if condition is met
                        if (aligned_counts[1] > min_reads) {
                                write.table(sample_counts, file = file.path(output_passed_count_dir, basename(f)), 
                                           sep = "\t", row.names = FALSE, quote = FALSE)
                        }
                }
        }
}

merge_counts <- function(output_passed_count_dir, output_merged_count_dir, input_library_count_file){

        # Listing files in the Passed Sample Counts directory
        dir_list <- list.files(path = file.path(output_passed_count_dir), pattern = "*.tsv", full.names = TRUE)

        # Looping through the files
        cores <-  detectCores()
        registerDoMC(cores)
        bh <- foreach (i = 1:length(dir_list)) %dopar% {
                f <- dir_list[i]
                if (!grepl("^\\.", basename(f))) { # Excluding hidden files
                        name <- tools::file_path_sans_ext(basename(f))
                        cat("Passed sample:",name, "\n")

                        # Reading the input library counts
                        input_counts <- read.csv2(file.path(input_library_count_file), 
                                                   sep = "\t", header = TRUE, col.names = c("id", "input"))

                        # Reading the sample counts
                        sample_counts <- read.csv2(f, sep = "\t", header = TRUE, check.names = FALSE)

                        # Merging the data
                        merged <- merge(input_counts, sample_counts, by = "id", all = TRUE)

                        # Writing the merged data to a file
                        merge_file <- paste0(name, "_merged.tsv")
			cat("Writing merged file:", merge_file, "\n")
                        write.table(merged, file = file.path(output_merged_count_dir, merge_file), 
                                    sep = "\t", row.names = FALSE, quote = FALSE)
                }
        }
}

user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}

# Define input folders
output_count_dir <- params$output_count_dir

input_library_count_file <- params$input_library_count_file

type <- params$type

# Create folders
if(type != "ntc"){
        output_passed_count_dir <- params$output_passed_count_dir
        if(dir.exists(output_passed_count_dir)){
                cat("Count directory already exists, output maybe over written.: ", output_passed_count_dir, "\n")
                user_input("Press [enter] to continue or [ctrl+z] to quit.")
        }else{
                cat("Creating output count dir: ", output_passed_count_dir, "\n")
                dir.create(output_passed_count_dir, recursive = TRUE, showWarnings = FALSE)
        }
}


output_merged_count_dir <- params$output_merged_count_dir
if(dir.exists(output_merged_count_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_merged_count_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_merged_count_dir, "\n")
        dir.create(output_merged_count_dir, recursive = TRUE, showWarnings = FALSE)
}

# Filter reads
if(type == "ntc"){
        merge_counts(output_count_dir, output_merged_count_dir, input_library_count_file)
} else{
        filter_reads(output_count_dir, output_passed_count_dir)
        merge_counts(output_passed_count_dir, output_merged_count_dir, input_library_count_file)
}



