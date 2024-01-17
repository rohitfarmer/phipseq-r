# PhiP-Seq Pipeline in R

## YAML file

### Input library
```yaml
# Inputs
type: input_library # [input_library, sample, ntc]. Sample inlucde study samples, anchors, and mocks.
demult_dir: "data/AACFV5MHV/Input"
index: "index/PanCov_PanicPotV2_Controls/PanCoV_PPV2_Controls"

# Outputs
output_count_dir: "results/input-library/counts"
output_trimmed_dir: "results/input-library/trimmed-data"

# Figures
figures_dir: "figures/input_library"

# Submit scripts and logs
submit_scr_dir: "results/input-library/submit-scripts"
hpc_log_dir: "results/input-library/hpc-log"

## Leave email blank if you prefer not to receive emails from the HPC
email: rohit.farmer@nih.gov
```

### Sample library
For all the other types of samples excluding the input library.

```yaml
# Inputs
type: sample # [input_library, sample, ntc]. Sample inlucde study samples, anchors, and mocks.
demult_dir: "data/AACFV5MHV/EVD68_Sample"
index: "index/PanCov_PanicPotV2_Controls/PanCoV_PPV2_Controls"
input_library_count_file: "results/input-library-counts/input-counts-avg.tsv"

# Outputs
## Replace "sample" with anchor/mocks/ntc to segregate output in their respective folders

output_count_dir: "results/sample/counts"
output_passed_count_dir: "results/sample/counts-passed"
output_merged_count_dir: "results/sample/counts-merged"
output_trimmed_dir: "results/sample/trimmed-data"
output_generalized_poisson_p_vals_dir: "results/sample/generalized-poisson-p-vals"
output_generalized_poisson_scores_dir: "results/sample/generalized-poisson-scores"
output_normalized_counts_dir: "results/sample/normalized-counts"

## Figures
figures_dir: "figures/sample"

# Submit scripts and logs
submit_scr_dir: "results/sample/submit-scripts"
hpc_log_dir: "results/sample/hpc-log"

## Leave email blank if you prefer not to receive emails from the HPC
email: rohit.farmer@nih.gov 
```

## Scripts
**Since the input and output folder paths are relative to the project root folder, always execute the scripts from the project root folder.**

### To perform sequence trimming and the Bowtie alignment
There is only one dependency of this script `yaml` which will be installed if it's not already install in your environment. 

```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-00-trim-bowtie.R meta/file.yaml
```

> For the scripts below **claim an interactive node** with as many cores as possible. Both the scripts can perform multithreaded parallel computation.

#### To plot overall Bowtie alignment rate (optional)
Fetch Bowtie alignment rate from the HPC log files to plot a histogram and save the data as a TSV file. This script only make sense for sample data.
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-01-bowtie-ali-rate-plot.R meta/sample_library.yaml
```

#### To plot pre-normalized counts distribution per sample (optional)
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-02-counts-dist.R meta/sample_library.yaml
```

#### To average counts if there are more than one input libraries (case based)
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-03-input-counts-avg.R meta/input_library.yaml
```

### To perform sample filter and merge
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr02-00-filter-merge.R meta/sample_library.yaml
```

### To perform phipstat normalization
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr03-00-phipstat-normalize.R
```
