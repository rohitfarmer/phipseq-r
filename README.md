# PhiP-Seq Pipeline in R

## YAML file

### Input library
```yaml
type: input_library
demult_dir: "data/Demultiplexed_Data/Input_Library"
output_count_dir: "results/Input_Library_Counts"
output_trimmed_dir: "results/Trimmed_Data"
index: "index/VirScan_Dereplicated_bowtie_index/VIR3_dereplicated"
submit_scr_dir: "results/submit-scripts"
hpc_log_dir: "results/hpc-log"
email: rohit.farmer@nih.gov
```

### Sample library
```yaml
# Inputs
type: sample
demult_dir: "data/Demultiplexed_Data/Sample_Data"
index: "index/VirScan_Dereplicated_bowtie_index/VIR3_dereplicated"
input_library_count_file: "results/Input_Library_Counts/VS_B1_3_counts.tsv"

# Outputs
output_count_dir: "results/Sample_Counts"
output_passed_count_dir: "results/Passed_Sample_Counts"
output_merged_count_dir: "results/Sample_Counts_Merged"
output_trimmed_dir: "results/Trimmed_Data"
output_generalized_poisson_p_vals_dir: "results/Generalized_Poisson_P_Vals"
output_generalized_poisson_scores_dir: "results/Generalized_Poisson_Scores"
output_normalized_sample_counts_dir: "results/Normalized_Sample_Counts"

# Submit scripts and logs
submit_scr_dir: "results/submit-scripts"
hpc_log_dir: "results/hpc-log"
email: rohit.farmer@nih.gov 
```

## Scripts
**Since the input and output folder paths are relative to the project root folder, always execute the scripts from the project root folder.**

### To perform sequence trimming and the Bowtie alignment
There is only one dependency of this script `yaml` which will be installed if it's not already install in your environment. 

```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr01-trim-bowtie.R meta/input.yaml
```

### To perform sample filter and merge
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr02-filter-merge.R meta/sample_library.yaml
```

### To perform phipstat normalization
```bash
module load R/4.3.1
Rscript --vanilla phipseq-r/scr03-phipstat-normalize.R
```
