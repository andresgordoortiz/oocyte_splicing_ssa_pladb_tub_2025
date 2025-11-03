#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=60

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=5G

# job name
#SBATCH --job-name vast-compare


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
singularity_image="docker://andresgordoortiz/vast-tools:latest"

# Run vast-tools compare using Singularity
singularity exec \
    --bind $PWD:/shared \
    $singularity_image bash -c "
        vast-tools compare \
            /shared/notebooks/inclusion_tables/pladb_INCLUSION_LEVELS_FULL-mm10.tab \
            -a 2022_038_S10_L001_R1_001_merged 2022_039_S11_L001_R1_001_merged 2022_040_S12_L001_R1_001_merged \
            -b 2022_044_S16_L001_R1_001_merged 2022_045_S17_L001_R1_001_merged 2022_046_S18_L001_R1_001_merged \
            --min_dPSI 10 \
            --min_range 3
    "


###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds


