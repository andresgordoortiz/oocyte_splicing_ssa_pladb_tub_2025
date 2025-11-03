#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=90

# queue
#SBATCH --qos=shorter
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name vast-align
# job array directive
#SBATCH --array=0-8

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

#Define file list and select the file for the current array job
files=($PWD/data/processed/pladb/*_merged.fastq.gz)
file=${files[$SLURM_ARRAY_TASK_ID]}


basename=$(basename "$file" .fastq.gz)
mkdir -p $PWD/data/processed/pladb/vast_out


singularity_image="docker://andresgordoortiz/vast-tools:latest"
VASTDB_PATH=$1
# Run vast-tools align using Singularity
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $PWD/data/processed/pladb/vast_out:/vast_out \
    $singularity_image vast-tools align \
    "$file" \
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
