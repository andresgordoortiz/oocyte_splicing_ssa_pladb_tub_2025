#!/usr/bin/bash


##################
# slurm setting #
##################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err

# time limit in minutes
#SBATCH --time=2880

# queue
#SBATCH --qos=vlong

# memory (MB)
#SBATCH --mem=36G
#SBATCH --cpus-per-task=2

# job name
#SBATCH --job-name run_RMarkdown

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

URL3="https://vastdb.crg.eu/downloads/mm10/PROT_IMPACT-mm10-v3.tab.gz"
FILE3="$PWD/notebooks/PROT_IMPACT-mm10-v2.3.tab.gz"
UNZIPPED_FILE3="${FILE3%.gz}"

if [ ! -f "$UNZIPPED_FILE3" ]; then
    if [ ! -f "$FILE3" ]; then
        echo "$FILE3 not found. Downloading..."
        wget "$URL3" -O "$FILE3"
    else
        echo "$FILE3 already exists. Skipping download."
    fi
    echo "Unzipping $FILE3..."
    gunzip -c "$FILE3" > "$UNZIPPED_FILE3"
else
    echo "$UNZIPPED_FILE3 already exists. Skipping download and unzip."
fi
# Set the working directory inside the container to /workspace
singularity run --bind "$(pwd)/notebooks:/shared" \
  docker://andresgordoortiz/splicing_analysis_r_crg:v1.5 \
  bash -c "cd /; Rscript -e \"rmarkdown::render('/shared/oocyte_tub_ssa_pladbnew.Rmd')\""


# end message #
###############
end_epoch=`date +%s`
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname) after $((end_epoch-start_epoch)) seconds
