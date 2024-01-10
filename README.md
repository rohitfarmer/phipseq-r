# PhiP-Seq Pipeline in R

## YAML file

### Input library
```yaml
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
demult_dir: "data/Demultiplexed_Data/Sample_Data"
output_count_dir: "results/Sample_Counts"
output_trimmed_dir: "results/Trimmed_Data"
index: "index/VirScan_Dereplicated_bowtie_index/VIR3_dereplicated"
submit_scr_dir: "results/submit-scripts"
hpc_log_dir: "results/hpc-log"
email: rohit.farmer@nih.gov 
```

## Scripts
**Since the input and output folder paths are relative to the project root folder, always execute the scripts from the project root folder.**

### To perform sequence trimming and the Bowtie alignment
There is only one dependency of this script `yaml` which will be installed if it's not already install in your environment. 

```bash
module load R
Rscript --vanilla phipseq-r/scr01-trim-bowtie.R meta/input.yaml
```
