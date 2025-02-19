#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')
if (!require('edgeR')) BiocManager::install("edgeR"); library("edgeR")

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr04-00-run-edger.R sample_library.yaml \n")
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

# Functions

run_edger <- function(sample_input_dat){
        group  <- ifelse(colnames(sample_input_dat) == "Input", "Input", "Sample")
        if(sum(group == "Sample") > 1){
                cat("Running EdgeR using more than one replicate.\n")
                dge <- DGEList(counts = sample_input_dat, group = group)
                dge <- calcNormFactors(dge)
                dge <- estimateDisp(dge)
                design <- model.matrix(~ group)
                fit <- glmFit(dge, design)
                lrt <- glmLRT(fit, coef = 2)  # coef=2 specifies the treatment vs control comparison
                results <- topTags(lrt, n = Inf, sort.by="none")
        }else if(sum(group == "Sample") == 1){
                cat("Running EdgeR with a single replicate.\n")
                dge <- DGEList(counts = sample_input_dat, group = group)
                dge <- calcNormFactors(dge)
                dge$common.dispersion <- 0.1  
                #dge <- estimateDisp(dge)
                design <- model.matrix(~ group)
                fit <- glmFit(dge, design)
                lrt <- glmLRT(fit, coef = 2)  # coef=2 specifies the treatment vs control comparison
                results <- topTags(lrt, n = Inf, sort.by="none")
        }

        return(results$table)
}

# Create output directories
output_edger <- params$output_edger_dir
if(dir.exists(output_edger)){
        cat("EdgeR results directory already exists, output maybe over written.: ", output_edger, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating EdgeR results dir: ", output_edger, "\n")
        dir.create(output_edger, recursive = TRUE, showWarnings = FALSE)
}

# Define input folders
sample_info <- read.table(params$sample_info_file, header = TRUE,
                          sep = "\t", check.names = FALSE)

output_passed_count_dir <- params$output_passed_count_dir

input_library_count_file <- params$input_library_count_file

input_count <- read.table(input_library_count_file, header = TRUE,
                        sep = "\t", check.names = FALSE)

sample_info_split <- split(sample_info, interaction(sample_info$subject_id, sample_info$visit))

for(i in 1:length(sample_info_split)){
        subject_id <- unique(sample_info_split[[i]]$subject_id)
        visit <- unique(sample_info_split[[i]]$visit)
 
        cat(i, "Subject:", subject_id, "Visit:", visit, "\n")

        dat <- sample_info_split[[i]]
        sample_input_dat <- data.frame()
        for(j in 1:nrow(dat)){
                passed_count_file <- file.path(output_passed_count_dir, paste0(dat[j, "sample_id"], "_counts.tsv"))
                if(file.exists(passed_count_file)){
                        temp_dat <- read.table(passed_count_file, header = TRUE, sep = "\t",
                                check.names = FALSE)
                        if(j == 1 |  nrow(sample_input_dat) == 0){
                                sample_input_dat <- temp_dat
                        }else{
                                sample_input_dat <- merge(sample_input_dat, temp_dat, by = "id")
                        }
                }
        }

        if(ncol(sample_input_dat >= 1)){
                sample_input_dat <- merge(sample_input_dat, input_count, by = "id") 
                rownames(sample_input_dat) <- sample_input_dat$id
                sample_input_dat <- as.matrix(sample_input_dat[,2:ncol(sample_input_dat)])
                
                edger_res <- run_edger(sample_input_dat)
                edger_res <- cbind("id" = rownames(edger_res), data.frame(edger_res))

                output_file_name <- file.path(output_edger, paste0(subject_id, "_", visit, "_edger.tsv"))
                write.table(edger_res, file = output_file_name, 
                            sep = "\t", row.names = FALSE, quote = FALSE)
        }
}



