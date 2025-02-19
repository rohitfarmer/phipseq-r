#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr04-01-merge-edger-scores.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "evd68-full", "sample_library.yaml")) # For testing keep it commented in production

user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}



# Create output directories
output_edger_scores <- params$output_edger_scores_dir
if(dir.exists(output_edger_scores)){
        cat("EdgeR merged scores directory already exists, output maybe over written.: ", output_edger_scores, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating EdgeR merged scores dir: ", output_edger_scores, "\n")
        dir.create(output_edger_scores, recursive = TRUE, showWarnings = FALSE)
}

# Define input folders
output_edger <- params$output_edger_dir
sample_info <- read.table(params$sample_info_file, header = TRUE,
                          sep = "\t", check.names = FALSE)

files_to_merge <- unique(paste0(sample_info$subject_id, "_", sample_info$visit))

df_out <- data.frame()
for(i in 1:length(files_to_merge)){
        cat("Reading:",i, files_to_merge[i], "\n")
        file_name <- file.path(output_edger, paste0(files_to_merge[i], "_edger.tsv")) 
        edger_res <- read.table(file_name, sep = "\t", header = TRUE, check.names = FALSE)
        edger_res <- edger_res[, c("id", "FDR")]
        colnames(edger_res) <- c("id", "PValue")
        edger_res$PValue <- -log10(edger_res$PValue)
        colnames(edger_res) <- c("id", files_to_merge[i])
        if(ncol(df_out) == 0){
                df_out <- edger_res
        }else if(ncol(df_out) > 0){
                df_out <- merge(df_out, edger_res, by = "id")
        }
}

write.table(df_out, file = file.path(output_edger_scores, "edger_scores_fdr.tsv"), 
                    sep = "\t", row.names = FALSE, quote = FALSE)


