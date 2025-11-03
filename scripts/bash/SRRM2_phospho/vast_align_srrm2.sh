#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=240
# queue
#SBATCH --qos=short
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=80G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name vast-align-paired
# job array directive - 9 pairs
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


###############
# run command #
###############

# Change to processed directory to get list of sample base names
processed_dir="$PWD/data/processed/SRRM2_phospho"
trimmed_dir="$processed_dir/trimmed"

# Get a list of all unique sample base names (without _1 or _2 suffix)
cd $trimmed_dir
mapfile -t base_names < <(ls *_val_1.fq.gz | sed 's/_1_val_1\.fq\.gz$//')

# Exit if array index is out of bounds
if [ $SLURM_ARRAY_TASK_ID -ge ${#base_names[@]} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of pairs (${#base_names[@]})"
    exit 1
fi

# Get base name corresponding to this array task
current_base=${base_names[$SLURM_ARRAY_TASK_ID]}

# Set file names for the pair of trimmed files
file1="${current_base}_1_val_1.fq.gz"
file2="${current_base}_2_val_2.fq.gz"

# Verify both files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: One or both paired files are missing ($file1, $file2)"
    exit 1
fi

# Create output directory
mkdir -p $processed_dir/vast_out

echo "Processing paired-end files: $file1 and $file2"

# Full paths to files
file1_path="$trimmed_dir/$file1"
file2_path="$trimmed_dir/$file2"

singularity_image="docker://andresgordoortiz/vast-tools:latest"
VASTDB_PATH=$1

# Run vast-tools align using Singularity in paired-end mode
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $processed_dir/vast_out:/vast_out \
    --bind $trimmed_dir:/trimmed \
    $singularity_image vast-tools align \
    "/trimmed/$file1" \
    "/trimmed/$file2" \
    -sp hg19 \
    -o /vast_out \
    --IR_version 2 \
    -c 8 \
    -n "$current_base"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds - task ID: $SLURM_ARRAY_TASK_ID