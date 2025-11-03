# FMN2/SPIRE Alternative Splicing Analysis Pipeline

**Author**: Andrés Gordo Ortiz
**Institution**: Centre for Genomic Regulation (CRG)
**Supervisors**: Adel Al Jord, Manuel Irimia
**Project Type**: Master Thesis

## Project Overview

This repository contains a comprehensive Nextflow pipeline for RNA-Seq data processing and alternative splicing analysis for the ***Al Jord Lab @CRG***. The pipeline enables the analysis of splicing events in mouse oocytes across multiple experimental conditions, including:

- *Pladienolide B* treatment at different doses (1mM and 10mM)
- *Double knock-out of FMN2 F-actin nucleating factor*
- *Double knock-out of the FMN2-interacting Spire protein*
- Additional datasets (SSA and TUB)

The pipeline streamlines data exploration, processing, and plotting for differential analysis of mRNA splicing events.

## Repository Structure

- **config/**: Configuration files for Nextflow, Docker, and cluster settings
- **data/**: Raw, processed and metadata files
  - **sample_sheets/**: CSV files defining samples for each experiment
- **scripts/**:
  - **R/notebooks/**: R Markdown notebooks for downstream analysis
- **main.nf**: The main Nextflow pipeline script
- **nextflow.config**: Configuration for the Nextflow pipeline
- **docs/**: Documentation files

## Nextflow Pipeline Features

Our unified Nextflow pipeline (`main.nf`) provides:

- **Reproducibility**: Containerized execution with Docker/Singularity
- **Scalability**: Works efficiently with both small and large datasets
- **Flexibility**: Configurable for various experimental designs via CSV sample sheets
- **Modularity**: Optional preprocessing steps (FastQC, trimming) that can be included/excluded
- **Integration**: Complete workflow from raw reads to alternative splicing analysis with VAST-tools

## Processing Workflow

The pipeline executes the following steps:

1. **Quality Control**: FastQC analysis of raw reads
2. **Read Trimming** (optional): Trim Galore for adapter and quality trimming
3. **VAST-tools Alignment**: Mapping to the mouse genome (mm10) and splicing event quantification
4. **Result Combination**: Merging of individual sample results
5. **Reporting**:
   - MultiQC report for QC metrics
   - R Markdown analysis report for splicing patterns

## Installation and Setup

### Requirements

- Nextflow (21.10.0 or later)
- Docker or Singularity (for containerized execution)
- VAST-tools database for your species of interest

### Running the Pipeline

1. Clone the repository:
   ```bash
   git clone https://github.com/andresgordoortiz/oocyte_fmn2_spire_pladb_crg.git
   cd oocyte_fmn2_spire_pladb_crg
   ```

2. Create a sample sheet CSV file (see examples in `data/sample_sheets/`)
   ```csv
   sample,fastq_1,fastq_2,type,group
   sample1,path/to/sample1_R1.fastq.gz,path/to/sample1_R2.fastq.gz,paired,control
   sample2,path/to/sample2_R1.fastq.gz,path/to/sample2_R2.fastq.gz,paired,treatment
   ```

3. Run the pipeline:
   ```bash
   nextflow run main.nf \
     --sample_csv data/sample_sheets/your_experiment.csv \
     --vastdb_path /path/to/vastdb \
     --data_dir /path/to/fastq_files \
     --outdir results/your_experiment \
     --species mm10
   ```

### Pipeline Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--sample_csv` | CSV file with sample information | (required) |
| `--vastdb_path` | Path to VAST-tools database | (required) |
| `--data_dir` | Directory containing FASTQ files | (required) |
| `--outdir` | Output directory | nextflow_results |
| `--species` | Species for alignment (e.g., mm10) | mm10 |
| `--skip_fastqc` | Skip FastQC step | false |
| `--skip_trimming` | Skip trimming step | false |
| `--skip_rmarkdown` | Skip R Markdown report | false |
| `--rmd_file` | Path to custom R Markdown | scripts/R/notebooks/Oocyte_fmndko_spireko_complete.Rmd |

### Running on a Cluster

For SLURM execution:

```bash
nextflow run main.nf \
  -profile slurm \
  --sample_csv data/sample_sheets/your_experiment.csv \
  --vastdb_path /path/to/vastdb \
  --data_dir /path/to/fastq_files
```

## Analyzing Different Datasets

The pipeline can be applied to all experimental conditions by using different sample sheets:

```bash
# For FMN2 knockout analysis
nextflow run main.nf --sample_csv data/sample_sheets/fmndko_samples.csv [other parameters]

# For Spire knockout analysis
nextflow run main.nf --sample_csv data/sample_sheets/spire_samples.csv [other parameters]

# For Pladienolide B treatment analysis
nextflow run main.nf --sample_csv data/sample_sheets/pladb_samples.csv [other parameters]
```

## Docker Containers

The pipeline utilizes the following Docker containers:
1. **VAST-tools**: `andresgordoortiz/vast-tools:latest`
2. **R/Markdown Processing**: `andresgordoortiz/splicing_analysis_r_crg:v1.5`
3. **FastQC**: `quay.io/biocontainers/fastqc:0.11.9--0`
4. **Trim Galore**: `https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0`
5. **MultiQC**: `multiqc/multiqc:latest`

## References
<a id="1">[1]</a>
Philip Ewels, Måns Magnusson, Sverker Lundin, Max Käller, MultiQC: summarize analysis results for multiple tools and samples in a single report, Bioinformatics, Volume 32, Issue 19, October 2016, Pages 3047–3048, https://doi.org/10.1093/bioinformatics/btw354

<a id="2">[2]</a>
Gohr, A., Mantica, F., Hermoso-Pulido, A., Tapial, J., Márquez, Y., Irimia, M. (2022). Computational Analysis of Alternative Splicing Using VAST-TOOLS and the VastDB Framework. In: Scheiffele, P., Mauger, O. (eds) Alternative Splicing. Methods in Molecular Biology, vol 2537. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-2521-7_7

<a id="3">[3]</a>
Ascensão-Ferreira, M., Martins-Silva, R., Saraiva-Agostinho, N. & Barbosa-Morais, N. L. betAS: intuitive analysis and visualization of differential alternative splicing using beta distributions. RNA 30, 337 (2024).