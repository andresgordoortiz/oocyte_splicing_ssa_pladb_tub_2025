#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=120

# queue
#SBATCH --qos=shorter
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=15G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name kallisto-pseudoalign

# array job settings
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

echo "Starting Kallisto pseudoalignment job"

# Define paths
FASTQ_DIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/data/processed/pladb"
OUTDIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Kallisto"
INDEX_DIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Kallisto/index"

# Get the list of FASTQ files
FASTQ_FILES=($(ls "$FASTQ_DIR"/*.fastq.gz))
FASTQ=${FASTQ_FILES[$((SLURM_ARRAY_TASK_ID-1))]}
SAMPLE_NAME=$(basename "$FASTQ" .fastq.gz)
SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_NAME"

# Create output directory for the sample
mkdir -p "$SAMPLE_OUTDIR"

#######################
# Run Kallisto pseudo #
#######################
echo "Processing $SAMPLE_NAME..."

singularity exec --bind $PWD:$PWD docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1 kallisto quant -i "$INDEX_DIR/kallisto_index" \
                  -o "$SAMPLE_OUTDIR" \
                  -t 8 \
                  --single -l 200 -s 20 "$FASTQ"

echo "Sample $SAMPLE_NAME processed successfully."
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds