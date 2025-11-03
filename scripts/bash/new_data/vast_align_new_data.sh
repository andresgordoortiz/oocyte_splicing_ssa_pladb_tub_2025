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
#SBATCH --job-name vast-align-single
# job array directive
#SBATCH --array=0-14

#################
# start message #
#################
start_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)

##################################
# make bash behave more robustly #
##################################
set -e
set -u
set -o pipefail


###############
# run command #
###############

# Define file list and select single file for the current array job
file_list=($PWD/data/processed/new_data/trimmed/*.fq.gz)

# Exit if array index is out of bounds
if [ $SLURM_ARRAY_TASK_ID -ge ${#file_list[@]} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds number of files (${#file_list[@]})"
    exit 1
fi

current_file=${file_list[$SLURM_ARRAY_TASK_ID]}

# Extract base name without extension
basename=$(basename "$current_file" _merged_trimmed.fq.gz)
mkdir -p $PWD/data/processed/new_data/vast_out

echo "Processing single-end file: $current_file"

singularity_image="docker://andresgordoortiz/vast-tools:latest"
VASTDB_PATH=$1

# Run vast-tools align using Singularity in single-end mode
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $PWD/data/processed/new_data/vast_out:/vast_out \
    $singularity_image vast-tools align \
    "$current_file" \
    -sp mm10 \
    -o /vast_out \
    --IR_version 2 \
    -c 8 \
    -n "$basename"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
