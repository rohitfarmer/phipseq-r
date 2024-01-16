#!/usr/bin/env Rscript


if (!require('yaml')) install.packages('yaml'); library('yaml')
if (!require('doMC')) install.packages('doMC'); library('doMC')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr01-02-counts-dist.R file.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "input_library.yaml")) # For testing keep it commented in production

user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}

# Create figures directory
figures_dir <- params$figures_dir
figures_dir <- file.path(figures_dir, "pre-normalized-counts-histograms")
if(dir.exists(figures_dir)){
        cat("Count directory already exists, output maybe over written.: ", figures_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", figures_dir, "\n")
        dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}


# Read count files
count_files <- list.files(file.path(params$output_count_dir), pattern = ".tsv",
                          full.names = TRUE)

type <- params$type

cores <-  detectCores()
registerDoMC(cores)
bh <- foreach(i = 1:length(count_files)) %dopar%{
        f <- count_files[i]
        df_counts <- read.csv2(f, sep = "\t", header = TRUE, check.names = FALSE)
        non_zero_counts <- df_counts[df_counts[2] != 0,][[2]]
        min_c <- min(non_zero_counts)
        max_c <- max(non_zero_counts)
        median_c <- median(non_zero_counts)
        mean_c <- mean(non_zero_counts)
        sd_c <- sd(non_zero_counts)
        range_c <- mean_c + sd_c
        plot_dat <- non_zero_counts[non_zero_counts <= range_c]
        
        # Plot the histogram
        plot_file <- file.path(figures_dir, paste0(tools::file_path_sans_ext(basename(f)), ".png"))
        cat("Generating:", i, ":", plot_file, "\n")
        png(filename = plot_file, width = 7 * 300, height = 5 * 300, res = 300)
        hist(plot_dat,                    
             main = paste0(tools::file_path_sans_ext(basename(f)), 
                           "\nnon-zero-counts: ", length(non_zero_counts),"/",nrow(df_counts),
                           "\nmin: ", min_c, 
                           "; max: ", max_c,
                           "; mean: ", round(mean_c,2),
                           "; median: ", median_c,
                           "; sd: ", round(sd_c,2)),
             xlab = "Counts (within one SD)",
             ylab = "Frequency",
             col = "#009E73",
             border = "black")
        abline(v = c(mean_c, median_c), col = c('black', 'blue'), lwd = 2, lty = 'dashed')

        dev.off()
}

