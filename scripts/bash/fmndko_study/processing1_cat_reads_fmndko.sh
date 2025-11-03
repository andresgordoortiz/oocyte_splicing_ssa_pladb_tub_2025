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
mkdir -p $PWD/data/processed/fmndko

input_dir="$PWD/data/raw/fmndko"
output_dir="$PWD/data/processed/fmndko"

# List all fastq.gz files in the input directory
files=($(ls "$input_dir"/*.fastq.gz))

# Iterate over the files in triples
for ((i=0; i<${#files[@]}; i+=3)); do
  # Define the output file name based on the first file in the triple
  output_file="$output_dir/$(basename ${files[i]} .fastq.gz)_$(basename ${files[i+1]} .fastq.gz)_$(basename ${files[i+2]} .fastq.gz)_merged.fastq.gz"

  # Concatenate the triple of files
  cat "${files[i]}" "${files[i+1]}" "${files[i+2]}" > "$output_file"

  # Print a message indicating the files have been merged
  echo "Merged ${files[i]}, ${files[i+1]}, and ${files[i+2]} into $output_file"
done

###############
# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds

