#!/usr/bin/env bash

#SBATCH --no-requeue

#SBATCH --mem 6G

#SBATCH -p genoa64

#SBATCH --qos pipelines



# Configure bash

set -e          # exit immediately on error

set -u          # exit immidiately if using undefined variables

set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails



# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal

# to Nextflow so that all jobs launched by the pipeline will be cancelled too.

_term() {

        echo "Caught SIGTERM signal!"

        # Send SIGTERM to the entire process group to ensure Tower gets the signal

        kill -s SIGTERM -$pid 2>/dev/null || kill -s SIGTERM $pid

        wait $pid

}



trap _term TERM INT



# load Java module

module load Nextflow

module load Java



# Check if we're running in a SLURM environment

if [ -n "$SLURM_JOB_ID" ]; then

    echo "Running as SLURM job ID: $SLURM_JOB_ID"

    echo "SLURM_CONF is set to: $SLURM_CONF"

else

    echo "Warning: Not running as a SLURM job. Process distribution may be limited."

fi



# Ensure SLURM_CONF is available to child processes

export SLURM_CONF=${SLURM_CONF:-/etc/slurm/slurm.conf}

echo "Using SLURM configuration: $SLURM_CONF"



# Make absolutely sure Nextflow sees this is a SLURM environment

export NXF_EXECUTOR=slurm

export NXF_CLUSTER_SEED=531684



# Ensure all relevant SLURM environment variables are preserved and passed to Nextflow

export SLURM_EXPORT_ENV=ALL

# Preserve additional SLURM variables that Tower might need

export SLURM_JOB_ID SLURM_JOB_NAME SLURM_CLUSTER_NAME



# limit the RAM that can be used by nextflow

export NXF_JVM_ARGS="-Xms2g -Xmx5g -Dexecutor.name=slurm"



# Debug output

echo "Environment variables for Nextflow:"

env | grep -E 'SLURM|NXF'



# Check if any arguments are provided

if [ $# -eq 0 ]; then

    echo "ERROR: No workflow file specified."

    echo "Usage: sbatch submit_nf.sh path/to/workflow.nf [additional parameters]"

    exit 1

fi



# Extract the workflow file (should be first argument)

WORKFLOW_FILE="$1"

shift



# Check if the workflow file exists or is a valid Nextflow workflow name

if [[ ! -f "$WORKFLOW_FILE" && ! "$WORKFLOW_FILE" =~ ^[a-zA-Z0-9_-]+$ ]]; then

    echo "WARNING: The specified workflow file '$WORKFLOW_FILE' does not exist as a file."

    echo "Nextflow will attempt to resolve it as a named workflow or URL."

fi



echo "Workflow file/name: $WORKFLOW_FILE"

echo "Additional parameters: $@"



# Always add executor override but REMOVE the -profile crg parameter

echo "Running Nextflow with executor=slurm"

CMD="nextflow run -ansi-log false $WORKFLOW_FILE $@"

echo "Executing: $CMD"



# Run Nextflow in foreground when using Tower to ensure proper signal handling

if [[ "$*" == *"-with-tower"* ]]; then

    echo "Tower detected - running in foreground for proper signal handling"

    eval "$CMD"

    exit_code=$?

else

    eval "$CMD" & pid=$!

    # Wait for the pipeline to finish

    echo "Waiting for Nextflow process ${pid}"

    wait $pid

    exit_code=$?

fi



# Print executor information from the trace file if available

if [ -f "trace.txt" ]; then

    echo "Executor information from trace file:"

    head -n 1 trace.txt

    grep -m 1 "executor" trace.txt || echo "No executor info found in trace file"

fi



# Return the actual exit status

exit $exit_code
