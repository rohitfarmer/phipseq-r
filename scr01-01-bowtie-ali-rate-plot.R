#!/usr/bin/env Rscript

# Fetch the overall Bowtie alignment rate from the HPC log files and plot the distribution.

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr01-01-bowtie-ali-rate-plot.R file.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "sample_library.yaml")) # For testing keep it commented in production

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
if(dir.exists(figures_dir)){
        cat("Count directory already exists, output maybe over written.: ", figures_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", figures_dir, "\n")
        dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

# Rread log files and fetch alignment rate percentage
hpc_log_dir <- file.path(params$hpc_log_dir)
log_files <- list.files(hpc_log_dir)
log_files <- sort(log_files)

df_out <- data.frame()
for(i in 1:length(log_files)){
        f <- log_files[i]
        sample_id <- sub("\\.o.*$", "", sub("jid_", "", f))
        bowtie_log <- readLines(file.path(hpc_log_dir, f))
        ali_rate <- strsplit(bowtie_log[15], " ")
        ali_rate <- as.numeric(sub("%", "", ali_rate[[1]][1]))
        df_out <- rbind(df_out,
                        data.frame("sample_id" = sample_id, 
                                   "overall_ali_rate" = ali_rate))
}
df_out <- na.omit(df_out)

# Plot the histogram
plot_file <- file.path(figures_dir, "overall-bowtie-ali-rate.png")
png(filename = plot_file, width = 7 * 300, height = 5 * 300, res = 300)
hist(df_out[["overall_ali_rate"]],
     main = "Overall Bowtie Alignment Rate",
     xlab = "Alignment Percentage",
     ylab = "Frequency",
     col = "blue",
     border = "black")
dev.off()

