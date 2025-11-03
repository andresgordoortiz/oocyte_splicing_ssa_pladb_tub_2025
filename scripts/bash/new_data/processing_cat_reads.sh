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
#SBATCH --mem=6G

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
mkdir -p $PWD/data/processed/new_data

input_dir="$PWD/data/raw/new_data"
output_dir="$PWD/data/processed/new_data"

# Get list of all "-a.fastq.gz" files
files_a=($(ls "$input_dir"/*-a.fastq.gz))

# Iterate over all "-a" files and find matching "-b" files
for file_a in "${files_a[@]}"; do
    # Construct the matching "-b" filename
    base_name=${file_a%-a.fastq.gz}
    file_b="${base_name}-b.fastq.gz"

    # Check if matching "-b" file exists
    if [ -f "$file_b" ]; then
        # Construct output filename (remove the "-a" suffix)
        sample_name=$(basename "$base_name" | tr '-' '_')
        output_file="$output_dir/${sample_name}_merged.fastq.gz"

        # Concatenate the pair of files
        cat "$file_a" "$file_b" > "$output_file"

        # Print a message indicating the files have been merged
        echo "Merged $file_a and $file_b into $output_file"
    else
        echo "Warning: No matching -b file found for $file_a"
    fi
done

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds