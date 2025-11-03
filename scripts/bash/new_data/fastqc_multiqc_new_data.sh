#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=60

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name=trim_fastqc_single

# Array job - process 15 files, max 5 concurrent jobs
#SBATCH --array=0-14

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname) - task ID: $SLURM_ARRAY_TASK_ID

##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail

######################
# create directories #
######################
mkdir -p $PWD/data/processed/new_data/trimmed
mkdir -p $PWD/data/processed/new_data/fastqc

#############################
# get files based on array ID #
#############################
# Store the original working directory
original_dir="$PWD"

# Change to processed directory - looking directly at processed files
processed_dir="$PWD/data/processed/new_data"
cd $processed_dir

# Get a list of all fastq.gz files (single-end)
mapfile -t fastq_files < <(ls *.fastq.gz)

# Exit if array index is out of bounds
if [ $SLURM_ARRAY_TASK_ID -ge ${#fastq_files[@]} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of files (${#fastq_files[@]})"
    exit 1
fi

# Get file corresponding to this array task
current_file=${fastq_files[$SLURM_ARRAY_TASK_ID]}

echo "Processing single-end file: $current_file"

# Run trim_galore in single-end mode using absolute paths
echo "Trimming file..."
singularity exec --bind "$processed_dir:/data/processed/new_data" \
    https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0 \
    trim_galore "/data/processed/new_data/$current_file" \
    --fastqc -j 8 -o /data/processed/new_data/trimmed -q 20 \
    --fastqc_args "-t 8 --outdir /data/processed/new_data/fastqc"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds - task ID: $SLURM_ARRAY_TASK_ID