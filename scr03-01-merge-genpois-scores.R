#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr04-01-merge-enrichment-scores.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "sample_library_mystat.yaml")) # For testing keep it commented in production


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
output_generalized_poisson_p_vals_dir <- params$output_generalized_poisson_p_vals_dir


# Create folders
output_generalized_poisson_scores_dir  <- params$output_generalized_poisson_scores_dir 

if(dir.exists(output_generalized_poisson_scores_dir)){
        cat("Output directory already exists, output maybe over written.: ", output_generalized_poisson_scores_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output dir: ", output_generalized_poisson_scores_dir, "\n")
        dir.create(output_generalized_poisson_scores_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the input_dat from the input file
files <- list.files(path = output_generalized_poisson_p_vals_dir, full.names = TRUE, pattern = "*.tsv")

# Read files and merge the scores in a single dataframe
df_scores_out <- data.frame()
for(file_i in 1:length(files)){
        input <- files[file_i]
        cat("Processing:",file_i,  basename(input), "\n")
        input_dat <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
        input_dat[2] <- round(input_dat[2],2)

        if(file_i == 1){
                df_scores_out <- input_dat
        } else {
                df_scores_out <- merge(df_scores_out, input_dat, by = "id", all = TRUE)
                cat("New dimensions:", dim(df_scores_out), "\n")
        }
}

output_file_name <- file.path(output_generalized_poisson_scores_dir, "general_mlxp.tsv")
write.table(df_scores_out, file = output_file_name, 
                    sep = "\t", row.names = FALSE, quote = FALSE)

