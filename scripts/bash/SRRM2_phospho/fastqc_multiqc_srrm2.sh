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
#SBATCH --job-name=trim_fastqc_paired

# Array job - process 9 pairs (18 files), each array job handles one pair
#SBATCH --array=0-8

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
mkdir -p $PWD/data/processed/SRRM2_phospho/trimmed
mkdir -p $PWD/data/processed/SRRM2_phospho/fastqc

#############################
# get files based on array ID #
#############################
# Store the original working directory
original_dir="$PWD"

# Change to processed directory - looking directly at processed files
raw_dir="$PWD/data/raw/SRRM2_phospho"
processed_dir="$PWD/data/processed/SRRM2_phospho"
cd $raw_dir

# Get a list of all unique sample base names (without _1 or _2 suffix)
mapfile -t base_names < <(ls *_1.fastq.gz | sed 's/_1\.fastq\.gz$//')

# Exit if array index is out of bounds
if [ $SLURM_ARRAY_TASK_ID -ge ${#base_names[@]} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of pairs (${#base_names[@]})"
    exit 1
fi

# Get base name corresponding to this array task
current_base=${base_names[$SLURM_ARRAY_TASK_ID]}

# Set file names for the pair
file1="${current_base}_1.fastq.gz"
file2="${current_base}_2.fastq.gz"

# Verify both files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: One or both paired files are missing ($file1, $file2)"
    exit 1
fi

echo "Processing paired-end files: $file1 and $file2"

# Run trim_galore in paired-end mode using absolute paths
echo "Trimming files..."
singularity exec --bind "$raw_dir:/data/raw/SRRM2_phospho,$processed_dir:/data/processed/SRRM2_phospho" \
    https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0 \
    trim_galore --paired \
    "/data/raw/SRRM2_phospho/$file1" \
    "/data/raw/SRRM2_phospho/$file2" \
    --fastqc -j 8 -o /data/processed/SRRM2_phospho/trimmed -q 20 \
    --fastqc_args "-t 8 --outdir /data/processed/SRRM2_phospho/fastqc"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds - task ID: $SLURM_ARRAY_TASK_ID