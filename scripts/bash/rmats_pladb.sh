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
#SBATCH --requeue

# memory (MB)
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# job name
#SBATCH --job-name rmats


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



singularity_image="docker://xinglab/rmats:v4.3.0"

# Run vast-tools align using Singularity
singularity exec --bind $PWD/STAR/index:/index \
    --bind $PWD/data/processed/pladb:/pladb \
    $singularity_image python /rmats/rmats.py \
    --s1 $PWD/data/processed/pladb/control.txt \
    --s2 $PWD/data/processed/pladb/pladb10mm.txt \
    --gtf $PWD/Mus_musculus.GRCm39.113.gtf \
    --bi $PWD/STAR/index \
    -t paired \
    --readLength 60 \
    --nthread 8 \
    --od $PWD/rmats_out \
    --tmp $PWD/tmp


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
