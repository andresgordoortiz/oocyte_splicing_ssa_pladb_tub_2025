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
#SBATCH --cpus-per-task=4

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

PROJECT_NAME=$1
DOWNLOADED_PATH=$2

if [ "$DOWNLOADED_PATH" != "0" ]; then
    FOLDER=$DOWNLOADED_PATH
else
    FOLDER=$PWD/data/raw/$PROJECT_NAME
fi

mkdir -p $FOLDER/fastqc
singularity exec --bind $FOLDER docker://biocontainers/fastqc:v0.11.9_cv8 fastqc -t 4 -o $FOLDER/fastqc $FOLDER/*.{fastq.gz,fq.gz}

################
# run multiqc  #
################
singularity exec --bind $FOLDER/fastqc docker://multiqc/multiqc:latest multiqc .

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds