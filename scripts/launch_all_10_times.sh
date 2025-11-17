#!/bin/bash

ROOT=$(cd "$(dirname "$0")/.." && pwd)


# This is the temporary folder that 'launch_all.sh' uses for it's output files
# (it MUST CORRISPOND TO "RESULTS_DIR" in 'launch_all.sh')
CURRENT_RUN_DIR="i_miei_risultati"
LAUNCH_SCRIPT="$ROOT/scripts/launch_all.sh"
FINAL_RESULTS_BASE_DIR="$ROOT/results"
EXPECTED_FILES=80
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
    file_count=0
    wait_time_files=0
    MAX_WAIT_FILE=300

    while true; do
        job_count=$(qstat -u $USER | wc -l)
        file_count=$(find "$CURRENT_RUN_DIR" -name "*.out" | wc -l)
        
        
        if [ "$job_count" -eq 0 ] && [ "$file_count" -ge "$EXPECTED_FILES" ]; then
            echo "Success!"
            break
        fi

        
        if [ "$job_count" -eq 0 ] && [ "$file_count" -lt "$EXPECTED_FILES" ]; then
            echo "Waiting for other files..."
            
            wait_time_files=$((wait_time_files + 30))
            
            # Se aspettiamo troppo (5 minuti a vuoto), forziamo l'uscita
            if [ "$wait_time_files" -ge "$MAX_WAIT_FILES" ]; then
                echo "TIMEOUT"
                break
            fi
        else
            wait_time_files=0
            echo "In queue: $job_count | File Ready: $file_count/$EXPECTED_FILES"
        fi

        sleep 30
    done

    echo "Batch $i completed!!"

    FINAL_RUN_NAME="run_${i}"
    
    mv "$CURRENT_RUN_DIR" "$FINAL_RESULTS_BASE_DIR/$FINAL_RUN_NAME"

    if [ "$file_count" -lt "$EXPECTED_FILES" ]; then
        echo "WARNING: There should be more files... Check manually :)"
    fi
done

echo "ALL COMPLETE!!"
