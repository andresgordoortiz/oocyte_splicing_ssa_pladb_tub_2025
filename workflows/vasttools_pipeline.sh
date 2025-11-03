#!/usr/bin/bash

##############################
# Full Processing Pipeline of Vast Tools CRG Adel Lab
# This script processes data by performing the following steps:
# 1. (optional) Download data from ENA
# 2. Align reads
# 3. Generate a multiQC report
# 5. Run vast combine
##############################

# SLURM output and error files
#SBATCH --output=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.out
#SBATCH --error=/users/aaljord/agordo/git/24CRG_ADEL_MANU_OOCYTE_SPLICING/logs/%x.%A_%a.err



# Input file with all samples to be downloaded
SAMPLES_FILE=$1
# Calculate the number of samples (number of rows)
NUM_SAMPLES=$(wc -l < "$SAMPLES_FILE")

PROJECT_NAME=$2

SPECIES=$3

# Ensure SPECIES is either "human" or "mouse"
if [[ "$SPECIES" != "hg38" && "$SPECIES" != "mm10" ]]; then
    echo "Error: SPECIES must be either 'hg38' or 'mm10'."
    exit 1
fi

echo "Species is set to $SPECIES"

VASTDB_PATH=$4
# Check if the user specified the read type, default to "single-end"
READ_TYPE=$5

# Ensure READ_TYPE is either "single-end" or "paired-end"
if [[ "$READ_TYPE" != "single-end" && "$READ_TYPE" != "paired-end" ]]; then
    echo "Error: READ_TYPE must be either 'single-end' or 'paired-end'."
    exit 1
fi

echo "Read type is set to $READ_TYPE"

# How to use: --skip-download to skip the download step, --downloaded-path to specify the path to the downloaded files
SKIP_DOWNLOAD=0
DOWNLOADED_PATH="0"

# Check if the last argument is --skip-download
if [[ "${!#}" == "--skip-download" ]]; then
    SKIP_DOWNLOAD=1
    DOWNLOADED_PATH=$SAMPLES_FILE
fi


# Debug: Print the parsed arguments
echo "SKIP_DOWNLOAD: $SKIP_DOWNLOAD"
echo "DOWNLOADED_PATH: $DOWNLOADED_PATH"

# Check if SKIP_DOWNLOAD is set and DOWNLOADED_PATH is provided
if [[ $SKIP_DOWNLOAD -eq 1 && -n "$DOWNLOADED_PATH" ]]; then
    echo "Skipping first job: Using pre-downloaded files at $DOWNLOADED_PATH..."
    # Submit second job
    echo "Submitting second job: Align reads..."
    jid2=$(sbatch $PWD/scripts/bash/step2_vastalign_vasttools.sh $NUM_SAMPLES $READ_TYPE $PROJECT_NAME $VASTDB_PATH $SPECIES $DOWNLOADED_PATH | tr -cd '[:digit:].')
    echo "...second job ID is $jid2"
else
    # Submit first job
    echo "Submitting first job: Download fastq..."
    jid1=$(sbatch $PWD/scripts/bash/step1_downloadfastq_vasttools.sh $SAMPLES_FILE $PROJECT_NAME | tr -cd '[:digit:].')
    echo "...first job ID is $jid1"

    # Submit second job dependent on the first
    echo "Submitting second job: Align reads..."
    jid2=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/step2_vastalign_vasttools.sh $NUM_SAMPLES $READ_TYPE $PROJECT_NAME $VASTDB_PATH $SPECIES | tr -cd '[:digit:].')
    echo "...second job ID is $jid2"
fi

# Third job - generate multiQC report (dependent on first job)
if [[ $SKIP_DOWNLOAD -eq 1 && -n "$DOWNLOADED_PATH" ]]; then
    echo "Submitting third job: Generate multiQC report..."
    jid3=$(sbatch $PWD/scripts/bash/step3_multiqc_vasttools.sh $PROJECT_NAME $DOWNLOADED_PATH | tr -cd '[:digit:].')
else
    echo "Submitting third job: Generate multiQC report..."
    jid3=$(sbatch --dependency=afterok:$jid1 $PWD/scripts/bash/step3_multiqc_vasttools.sh $PROJECT_NAME $DOWNLOADED_PATH | tr -cd '[:digit:].')
fi
echo "...third job ID is $jid3"

# Fourth job - run vast combine (dependent on third job)
echo "Submitting fifth job: Run vast combine..."
jid4=$(sbatch --dependency=afterok:$jid2 $PWD/scripts/bash/step4_vastcombine_vasttools.sh $PROJECT_NAME $VASTDB_PATH $SPECIES | tr -cd '[:digit:].')
echo "...fifth job ID is $jid4"


echo "All jobs submitted!"