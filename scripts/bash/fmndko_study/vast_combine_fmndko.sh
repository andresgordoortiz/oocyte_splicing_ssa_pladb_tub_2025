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


###############
# run command #
###############

VASTDB_PATH=$1

# Define Singularity image path
singularity_image="docker://andresgordoortiz/vast-tools:latest"

# Run vast-tools align using Singularity
singularity exec --bind $VASTDB_PATH:/usr/local/vast-tools/VASTDB \
    --bind $PWD/data/processed/fmndko:/fmndko \
    $singularity_image bash -c "vast-tools combine /fmndko/vast_out/to_combine -sp mm10 -o /fmndko/vast_out"

mv $PWD/data/processed/fmndko/vast_out/INCLUSION_LEVELS_FULL* $PWD/notebooks/inclusion_tables/fmndko_INCLUSION_LEVELS_FULL-mm10.tab
###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
