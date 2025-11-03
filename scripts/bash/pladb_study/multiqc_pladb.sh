#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=50

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
# job name
#SBATCH --job-name fastqc_multiqc


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

mkdir -p $PWD/data/processed/pladb/fastqc
singularity exec --bind $PWD/data/processed/pladb \
    docker://biocontainers/fastqc:v0.11.9_cv8 \
    fastqc -t 8 -o $PWD/data/processed/pladb/fastqc \
    $PWD/data/processed/pladb/*.fastq.gz

################
# run multiqc  #
################
singularity exec --bind $PWD/data/processed/pladb:/pladb \
    docker://multiqc/multiqc:latest \
    /bin/bash -c "cd /pladb && multiqc . -n pladb_multiqc_report.html"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds