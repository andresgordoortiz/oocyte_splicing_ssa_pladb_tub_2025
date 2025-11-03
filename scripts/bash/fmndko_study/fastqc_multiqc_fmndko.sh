#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A.err

# time limit in minutes
#SBATCH --time=60

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=5G
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



################
# run fastqc   #
################


mkdir -p $PWD/data/processed/fmndko/fastqc
singularity exec --bind $PWD/data/processed/fmndko \
    docker://biocontainers/fastqc:v0.11.9_cv8 \
    fastqc -t 8 -o $PWD/data/processed/fmndko/fastqc \
    $PWD/data/processed/fmndko/*.{fastq.gz,fq.gz}

################
# run multiqc  #
################
singularity exec --bind $PWD/data/processed/fmndko:/fmndko \
    docker://multiqc/multiqc:latest \
    /bin/bash -c "cd /fmndko && multiqc . -n fmndko_multiqc_report.html"

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds