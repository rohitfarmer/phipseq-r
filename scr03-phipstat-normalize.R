#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr02-filter-merge.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "sample_library.yaml")) # For testing keep it commented in production

# Functions
norm_round <- function(normalized_data_file){

        # Reading the tab-separated file
        normalized_data <- read.csv2(normalized_data_file, sep = "\t", header = TRUE, check.names = FALSE)

        # Rounding the data
        normalized_data[,2] <- as.numeric(normalized_data[,2])
        normalized_data[,3] <- as.numeric(normalized_data[,3])
        normalized_data_rounded <- cbind(normalized_data["id"], round(normalized_data[2:3]))

        # Modifying specific columns to remove '.0'
        #normalized_data_rounded[, 2] <- gsub("\\.0$", "", as.character(normalized_data_rounded[, 2]))
        #normalized_data_rounded[, 3] <- gsub("\\.0$", "", as.character(normalized_data_rounded[, 3]))

        # Return rounded data
        return(normalized_data_rounded)
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
output_merged_count_dir <- params$output_merged_count_dir


# Create folders
output_generalized_poisson_p_vals_dir <- params$output_generalized_poisson_p_vals_dir
if(dir.exists(output_generalized_poisson_p_vals_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_generalized_poisson_p_vals_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_generalized_poisson_p_vals_dir, "\n")
        dir.create(output_generalized_poisson_p_vals_dir, recursive = TRUE, showWarnings = FALSE)
}

output_generalized_poisson_scores_dir <- params$output_generalized_poisson_scores_dir
if(dir.exists(output_generalized_poisson_scores_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_generalized_poisson_scores_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_generalized_poisson_scores_dir, "\n")
        dir.create(output_generalized_poisson_scores_dir, recursive = TRUE, showWarnings = FALSE)
}

output_normalized_sample_counts_dir <- params$output_normalized_sample_counts_dir
if(dir.exists(output_normalized_sample_counts_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_normalized_sample_counts_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_normalized_sample_counts_dir, "\n")
        dir.create(output_normalized_sample_counts_dir, recursive = TRUE, showWarnings = FALSE)
}

# Processing files in Sample_Counts_Merged directory
files <- list.files(path = output_merged_count_dir, full.names = TRUE, pattern = "*.tsv")

for (i in 1:length(files)) {
        f <- files[i]
        f <- tools::file_path_sans_ext(basename(f))

        # Normalizing counts
        cat("Normalizing sample: ", f, "\n")
        norm_count_cmd <- paste("module load phipstat\n", "phip normalize-counts -i", 
                                shQuote(file.path(files[i])), "-o", 
                                shQuote(file.path(output_normalized_sample_counts_dir, paste0(f, ".tsv"))), 
                                "-m col-sum")
        system(norm_count_cmd)
        
        normalized_data_file <- file.path(output_normalized_sample_counts_dir, paste0(f, ".tsv"))
        normalized_data_rounded <- norm_round(normalized_data_file)
        write.table(normalized_data_rounded, file = normalized_data_file, 
                    sep = "\t", row.names = FALSE, quote = FALSE)

        # Compute p-values
        cat("Computing p-values for sample: ", f, "\n")
        phip_p_cmd <- paste0("module load phipstat\n", "phip compute-pvals -i ", 
                             shQuote(normalized_data_file), " -o ", 
                             shQuote(file.path(output_generalized_poisson_p_vals_dir, paste0(f, "_mlxp.tsv"))))
        system(phip_p_cmd)
}

# Merging columns
cat("Merging phipstat results ...\n")
merge_cmd <- paste0("module load phipstat\n", "phip merge-columns -i ", output_generalized_poisson_p_vals_dir, " -o ", 
                    file.path(output_generalized_poisson_scores_dir, "general_mlxp.tsv")," -p 1")
system(merge_cmd)



