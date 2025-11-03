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

# memory (MB)
#SBATCH --mem=100G
#SBATCH --cpus-per-task=6

# job name
#SBATCH --job-name rMATS_BAM


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
singularity_image="docker://xinglab/rmats:v4.3.0"


# Run rMATS using Singularity
singularity exec \
    --bind $PWD:/share \
    $singularity_image \
    bash -c "cd /rmats && python rmats.py \
    --b1 /share/data/processed/pladb/control_pladb_rmats.txt \
    --b2 /share/data/processed/pladb/pladb10mm_pladb_rmats.txt \
    --gtf /share/Mus_musculus.GRCm39.113.gtf \
    --readLength 68 \
    --nthread 6 \
    --od /share/rMATS_bam \
    --novelSS \
    --allow-clipping --task both \
    --tmp /share/tmp"
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
