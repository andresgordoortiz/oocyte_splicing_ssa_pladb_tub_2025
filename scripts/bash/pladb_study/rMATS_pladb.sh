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

# memory (MB)
#SBATCH --mem=50G
#SBATCH --cpus-per-task=6

# job name
#SBATCH --job-name rMATS


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

# Define Singularity image path
singularity_image="docker://mcfonsecalab/rmats:4.1.2"

# Run rMATS using Singularity
singularity exec \
    --bind $PWD/data/processed/pladb:/pladb \
    $singularity_image \
    python /rmats-turbo/rmats.py \
    --s1 /pladb/control_pladb_rmats.txt \
    --s2 /pladb/pladb10mm_pladb_rmats.txt \
    --gtf $PWD/Mus_musculus.GRCm39.113.gtf \
    --bi $PWD/STAR/index \
    --readLength 50 \
    --nthread 6 \
    --od $PWD/rMATS \
    --t single \
    --novelSS \
    --tmp $PWD/tmp

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
