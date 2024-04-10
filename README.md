# PhiP-Seq Pipeline in R

## YAML file

The current form of the pipeline utilizes a YAML file to pass on the parameters to the scripts. Below are the formats for two types of YAML files: first for the pre-immunoprecipitated phage library (a.k.a. input library) and second for post-immunoprecipitated phage libraries, which includes study sample (serum samples), NTCs (no template controls for the PCR steps), Anchor (control serum samples that we use on every plate to evaluate batch effects), and Mocks (phages run through the entire immunoprecipitation protocol without serum added). 

**Note:** All the data types require a separate YAML file, such as `input.yaml`, `sample.yaml`, `ntc.yaml`, `anchor.yaml`, and `mocks.yaml`. Please copy the contents from the sections below to create the data type-specific YAML files. 

### Pre-immunoprecipitation data type (input library)
```yaml
# Inputs
type: input_library # [input_library, sample, ntc]. Sample inlucde study samples, anchors, and mocks.
demult_dir: "data/AACFV5MHV/Input"  # The terminal folder in the path should lead to sub folders per sample. These subfolders should have two FASTQ files with forward and reverse reads.
index: "index/PanCov_PanicPotV2_Controls/PanCoV_PPV2_Controls"

# Outputs
output_trimmed_dir: "results/input-library/trimmed-data"
output_count_dir: "results/input-library/counts"

# Figures
figures_dir: "figures/input_library"

# Submit scripts and logs
submit_scr_dir: "results/input-library/submit-scripts"
hpc_log_dir: "results/input-library/hpc-log"

## Leave email blank if you prefer not to receive emails from the HPC
email: rohit.farmer@nih.gov
```

### Post-immunoprecipitation data type (sample, anchors, NTCs, and mocks)
```yaml
# Inputs
type: sample # [input_library, sample, ntc]. Sample inlucde study samples, anchors, and mocks.
demult_dir: "data/AACFV5MHV/EVD68_Sample"
index: "index/PanCov_PanicPotV2_Controls/PanCoV_PPV2_Controls"
input_library_count_file: "results/input-library-counts/input-counts-avg.tsv"

# Outputs
## Replace "sample" with anchor/mocks/ntc to segregate output in their respective folders

output_trimmed_dir: "results/sample/trimmed-data"
output_count_dir: "results/sample/counts"
output_passed_count_dir: "results/sample/counts-passed"
output_merged_count_dir: "results/sample/counts-merged"
output_normalized_counts_dir: "results/sample/normalized-counts"
output_generalized_poisson_p_vals_dir: "results/sample/generalized-poisson-p-vals"
output_generalized_poisson_scores_dir: "results/sample/generalized-poisson-scores"

## Figures
figures_dir: "figures/sample"

# Submit scripts and logs
submit_scr_dir: "results/sample/submit-scripts"
hpc_log_dir: "results/sample/hpc-log"

## Leave email blank if you prefer not to receive emails from the HPC
email: rohit.farmer@nih.gov 
```

## Scripts
**Note: Since the input and output folder paths are relative to the project root folder, always execute the scripts from the project root folder.**

### To perform sequence trimming and the Bowtie alignment

```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-00-trim-bowtie.R meta/file.yaml
```

**Note:** For the scripts below *claim an interactive node* with as many cores as possible to take advantage of multithreaded computing. 

#### 1. To plot overall Bowtie alignment rate (QC)
Fetch Bowtie alignment rate from the HPC log files to plot a histogram and save the data as a TSV file. This script only make sense for sample data.
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-01-bowtie-ali-rate-plot.R meta/sample_library.yaml
```

#### 1.1. To average counts if there are more than one input libraries (case based)
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-02-input-counts-avg.R meta/input_library.yaml
```

#### 1.2. To plot pre-normalized counts distribution per sample (QC)
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-03-counts-dist.R meta/sample_library.yaml
```

### 2. To perform sample filtering and merge
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr02-00-filter-merge.R meta/sample_library.yaml
```

### 3. To perform phipstat normalization and generalized Poisson modeling (deprecated) 
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/phipstat/scr03-00-phipstat-normalize.R sample_library.yaml
```


**OR**


### 2.1. To perform normaliztion on filtered and merged data
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr02-01-normalize.R sample_library.yaml
```

### 3. To conduct generalized Poisson modeling 
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr03-00-genpoise-scores.R sample_library.yaml
```

#### 3.1. To merge per sample generalized Poisson modeling scores into a single dataframe
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr03-01-merge-genpois-scores.R sample_library.yaml
```




