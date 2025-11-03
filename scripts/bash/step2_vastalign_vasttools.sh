#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=360

# queue
#SBATCH --qos=short

# memory (MB)
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name vast-align

#SBATCH --array=0-11

NUM_SAMPLES=$1
READ_TYPE=$2
PROJECT_NAME=$3
VASTDB_PATH=$4
SPECIES=$5
DOWNLOADED_PATH=$6


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

# Define file list and select the file for the current array job
if [ -z "$DOWNLOADED_PATH" ]; then
    files=($PWD/data/raw/$PROJECT_NAME/*q.gz)
else
    files=($PWD/$DOWNLOADED_PATH/*q.gz)
fi

file1=${files[$SLURM_ARRAY_TASK_ID]}

basename=$(basename "$file1" | sed -E 's/\.(fastq|fq)\.gz$//')
mkdir -p $PWD/data/processed/$PROJECT_NAME/vast_out

# Define Singularity image path
singularity_image="docker://andresgordoortiz/vast-tools:latest"

# Run vast-tools align using Singularity
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $PWD/data/processed/$PROJECT_NAME/vast_out:/vast_out \
    $singularity_image vast-tools align \
    "$file1" \
    -sp $SPECIES \
    -o /vast_out \
    --IR_version 2 \
    -c 8 \
    -n "$basename"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
