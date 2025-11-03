#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=20

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=5G

# job name
#SBATCH --job-name vast-combine


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

PROJECT_NAME=$1
VASTDB_PATH=$2
SPECIES=$3
###############
# run command #
###############
# Store current working directory
current_dir=$PWD
cd $PWD/data/processed/$PROJECT_NAME/vast_out/to_combine

# Define Singularity image path
singularity_image="docker://andresgordoortiz/vast-tools:latest"

# Run vast-tools align using Singularity
singularity exec $singularity_image vast-tools combine \
    --bind $VASTDB_PATH:/VASTDB \
    -sp $SPECIES \
    -o $current_dir/data/processed/$PROJECT_NAME/vast_out

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
