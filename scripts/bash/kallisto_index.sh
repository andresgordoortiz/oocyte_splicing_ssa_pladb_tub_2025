#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=300

# queue
#SBATCH --qos=short
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name kallisto-index

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


# Define paths
FASTA="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Mus_musculus.GRCm39.cdna.all.fa"
GTF="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Mus_musculus.GRCm39.113.gtf"
FASTQ_DIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/data/processed/pladb"
OUTDIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Kallisto"
INDEX_DIR="/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/Kallisto/index"

# Path to Singularity image\NSINGULARITY_IMAGE="docker://quay.io/biocontainers/kallisto:0.48.0--h87f3376_3"

# Create output directories
mkdir -p "$INDEX_DIR" "$OUTDIR"

###########################
# Build Kallisto index if needed #
###########################

singularity exec --bind $PWD:$PWD docker://quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1 kallisto index -i "$INDEX_DIR/kallisto_index" "$FASTA"

# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds