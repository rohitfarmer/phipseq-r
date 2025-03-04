#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr01-trim-bowtie.R file.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "evd68-sky", "anchor.yaml")) # For testing keep it commented in production

# Functions
exec_trim_bowtie <- function(type, d, files, read1, read2, output_count_dir, output_trimmed_dir, submit_scr_dir, hpc_log_dir, email = NULL){

        dir.create(file.path(output_trimmed_dir, d), recursive = TRUE, showWarnings = FALSE)
        count_file <- paste0(d, "_counts.tsv")
        
        submit_file <- file.path(submit_scr_dir, paste0(d, ".sh"))
        sink(submit_file)
        cat("#!/bin/bash\n")

        cat("# Job name\n")
        cat(paste0("#SBATCH --job-name=", d,"\n"))

        cat("# Job (logs and screen output) output location and file name\n")
        cat(paste0('#SBATCH --output=', file.path(hpc_log_dir, paste0(d,".txt")), '\n'))

        cat("# Number of tasks\n")
        cat("#SBATCH --ntasks=1\n")

        cat("# Number of CPUs per task\n")
        cat("#SBATCH --cpus-per-task=4\n")

        cat("# Memory per CPU\n")
        cat("#SBATCH --mem-per-cpu=1G\n")


        if(!is.null(email)){
                cat("# Send email when job is submitted and completed\n")
                cat(paste0("#SBATCH --mail-user=", email, "\n"))
                cat(paste0("#SBATCH --mail-type=FAIL", "\n"))
        }

        cat("\n")

        cat("cd $PWD\n")

        cat("container=/data/vrc_his/douek_lab/premise-bioinfo/containers/phipseq/phipstat-pipeline.sif\n")

        # Trim using fastx trimmer
        cat('singularity exec --bind "$PWD" $container gunzip -c', files[1],  '| singularity exec --bind "$PWD" $container fastx_trimmer -f 21 -o', 
            file.path(output_trimmed_dir, d, tools::file_path_sans_ext(read1)))
        cat("\n")
        cat('singularity exec --bind "$PWD" $container gunzip -c', files[2],  '| singularity exec --bind "$PWD" $container fastx_trimmer -f 28 -o', 
            file.path(output_trimmed_dir, d, tools::file_path_sans_ext(read2)))

        cat("\n\n")

        if(type == "input_library"){
                cat('echo -e "id\\tInput" >', file.path(output_count_dir, count_file), "\n\n")
        }else if(type == "sample"){
                cat(paste0('echo -e "id\\t', d, '" > ', file.path(output_count_dir, count_file)), "\n\n")
        }

        # Run bowtie2 and samtools
        cat('singularity exec --bind "$PWD" $container bowtie2 -p 4 -x', bowtie_index, "-1", 
                   file.path(output_trimmed_dir, d, tools::file_path_sans_ext(read1)), "-2", 
                   file.path(output_trimmed_dir, d, tools::file_path_sans_ext(read2)),
                   '| singularity exec --bind "$PWD" $container samtools sort -O BAM | singularity exec --bind "$PWD" $container samtools depth -aa -m 100000000 -',
                   "| awk 'BEGIN {OFS=\"\\t\"} {counts[$1] = ($3 < counts[$1]) ? counts[$1] : $3} END {for (c in counts) {print c, counts[c]}}'",    
                   "| sort -k 1",
                   ">>", file.path(output_count_dir, count_file), "\n")

        sink()

        cat("Submitting the job ... ")
        system(paste0("sbatch ", submit_file))
        cat("\n")
}

user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}

# Define input directories
demult_dir <- file.path(params$demult_dir)
if(dir.exists(demult_dir)){
        cat("Demultiplexed data directory: ", demult_dir,  "\n")
}else{
        cat("Demultiplexed data directory does not exist. Please check the path.\n")
}

bowtie_index <- file.path(params$index)
if(dir.exists(dirname(bowtie_index))){
        cat("Bowtie index directory: ", bowtie_index,  "\n")
}else{
        cat("Bowtie index directory does not exist. Please check the path.\n")
}

# Create output directories
output_count_dir <- params$output_count_dir
if(dir.exists(output_count_dir)){
        cat("Count directory already exists, output maybe over written.: ", output_count_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output count dir: ", output_count_dir, "\n")
        dir.create(output_count_dir, recursive = TRUE, showWarnings = FALSE)
}

output_trimmed_dir <- params$output_trimmed_dir
if(dir.exists(output_trimmed_dir)){
        cat("Trimmed data directory already exists, output maybe over written.: ", output_trimmed_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating trimmed data dir: ", output_trimmed_dir, "\n")
        dir.create(output_trimmed_dir, recursive = TRUE, showWarnings = FALSE)
}

submit_scr_dir <- params$submit_scr_dir
if(dir.exists(submit_scr_dir)){
        cat("Submit script directory already exists, output maybe over written.: ", submit_scr_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating submit script dir: ", submit_scr_dir, "\n")
        dir.create(submit_scr_dir, recursive = TRUE, showWarnings = FALSE)
}

hpc_log_dir <- params$hpc_log_dir
if(dir.exists(hpc_log_dir)){
        cat("HPC log directory already exists, output maybe over written.: ", hpc_log_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating HPC log dir: ", hpc_log_dir, "\n")
        dir.create(hpc_log_dir, recursive = TRUE, showWarnings = FALSE)
}

# Process each directory in the read library
input_dirs <- list.dirs(demult_dir, full.names=TRUE, recursive=FALSE)
for(i in 1:length(input_dirs)){
        input_dir <- input_dirs[i]
        cat("Generating HPC submission script for:", input_dir, "\n")
        d <- basename(input_dir)

        # List all files in the directory
        files <- list.files(path=input_dir, full.names=TRUE)
        #cat(files, "\n")

        # Assuming there are two reads
        read1 <- basename(files[1])
        cat("Read 1: ", read1, "\n")
        read2 <- basename(files[2])
        cat("Read 2: ", read2, "\n")

        exec_trim_bowtie(params$type, d, files, read1, read2, output_count_dir, output_trimmed_dir, submit_scr_dir, hpc_log_dir, params$email)
}






