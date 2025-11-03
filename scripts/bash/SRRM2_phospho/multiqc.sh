#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=10

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1

# job name
#SBATCH --job-name=multiqc


################
# run multiqc  #
################
singularity exec --bind $PWD/data/processed/SRRM2_phospho:/SRRM2_phospho \
    docker://multiqc/multiqc:latest \
    /bin/bash -c "cd /SRRM2_phospho && multiqc . -n SRRM2_phospho_multiqc_report.html"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds