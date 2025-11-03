#!/bin/bash


##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=30

# queue
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=4G

# job name
#SBATCH --job-name cat_reads


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
# Define the input and output directories
mkdir -p $PWD/data/processed/pladb

input_dir="$PWD/data/raw/pladb"
output_dir="$PWD/data/processed/pladb"

# List all fastq.gz files in the input directory
files=($(ls "$input_dir"/*.fastq.gz))

# Iterate over the files in pairs
for ((i=0; i<${#files[@]}; i+=2)); do
    # Define the output file name based on the first file in the pair
    output_file="$output_dir/$(basename ${files[i]} .fastq.gz | tr '-' '_')_merged.fastq.gz"

    # Concatenate the pair of files
    cat "${files[i]}" "${files[i+1]}" > "$output_file"

    # Print a message indicating the files have been merged
    echo "Merged ${files[i]} and ${files[i+1]} into $output_file"
done

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

