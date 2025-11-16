#!/bin/bash

ROOT=$(cd "$(dirname "$0")/.." && pwd)


# This is the temporary folder that 'launch_all.sh' uses for it's output files
# (it MUST CORRISPOND TO "RESULTS_DIR" in 'launch_all.sh')
CURRENT_RUN_DIR="i_miei_risultati"


LAUNCH_SCRIPT="$ROOT/scripts/launch_all.sh"

FINAL_RESULTS_BASE_DIR="$ROOT/results"
mkdir -p "$FINAL_RESULTS_BASE_DIR"


cd "$REPO_ROOT"
if [ $? -ne 0 ]; then
    echo "ERROR: Can't change Directory in: $REPO_ROOT"
    exit 1
fi


echo "Starting..."
echo "Final results will be into: $FINAL_RESULTS_BASE_DIR"

for i in $(seq 1 10); do
    echo ""
    echo "Starting Batch: $i"

    #Cleaning the temporary folder
    rm -rf "$CURRENT_RUN_DIR"
    
    "$LAUNCH_SCRIPT"
    
    echo "Done, waiting for completion..."


    job_count=$(qstat -u $USER | wc -l)

    while [ $job_count -gt 0 ]; do
        echo "Still waiting for: $job_count jobs to terminate"
        sleep 60
        job_count=$(qstat -u $USER | wc -l)
    done

    echo "Batch $i completed!!"

    FINAL_RUN_NAME="run_${i}"
    
    mv "$CURRENT_RUN_DIR" "$FINAL_RESULTS_BASE_DIR/$FINAL_RUN_NAME"
done

echo "ALL COMPLETE!!"
