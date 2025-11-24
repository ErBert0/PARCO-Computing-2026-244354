#!/bin/bash

#1. This finds the root directory of the project
ROOT=$(cd "$(dirname "$0")/.." && pwd)

#2.This is the absolute path for the PBS template
PBS_TEMPLATE="$ROOT/scripts/job_template.pbs"
RESULTS_DIR="$ROOT/results_timings"


# --- MAIN VARIABLES TO CHANGE---
WALLTIME="03:00:00"
MEMORY=20gb
MAX_JOBS=30
#------------------------------------------------------------------------


# --- CONFIGURATION --- (Change Here directorys)
#Here you HAVE to specify the absolute path of the matrix folder
MATRICI_DIR="/home/matteo.berte/Matrices"
MATRICI=(
    "$MATRICI_DIR/1_ASIC.mtx"
    "$MATRICI_DIR/2_kkt.mtx"
    "$MATRICI_DIR/3_mawi.mtx"
    "$MATRICI_DIR/4_nlpkk.mtx"
    "$MATRICI_DIR/5_ecology2.mtx"
)
THREADS=(4 8 16 32 64) # List of threads to test
SCHEDULERS=("static" "dynamic" "guided") #Schedulers


mkdir -p $RESULTS_DIR


check_queue_limit() {
    # Contiamo quante righe 'matteo.b' appaiono in qstat
    local job_count=$(qstat -u $USER |wc -l) 
    
    while [ $job_count -ge $MAX_JOBS ]; do
        echo "  Queue full: ($job_count / $MAX_JOBS)"
        sleep 60
        job_count=$(qstat -u $USER | wc -l)
    done
}

# --- SEQUENTIAL SECTION ---
echo "Sending Sequential Jobs"
for M in "${MATRICI[@]}"; do
    M_NAME=$(basename "$M" .mtx)
    
    # Name file format: "matrix_scheduler_threads"
    JOB_NAME="${M_NAME}_sequential_1"
    OUT_FILE="${RESULTS_DIR}/${JOB_NAME}.out"
    ERR_FILE="${RESULTS_DIR}/${JOB_NAME}.err"
    
    check_queue_limit

    echo "Sending: $JOB_NAME"
    
    qsub -N "$JOB_NAME" \
         -o "$OUT_FILE" \
         -e "$ERR_FILE" \
         -l select=1:ncpus=1:mem=$MEMORY \
         -l walltime=$WALLTIME \
         -v ROOT="$ROOT",MATRICE="$M",SCHEDULER="sequential",NCPUS=1 \
         "$PBS_TEMPLATE"
    sleep 1
done

# --- 2. PARALLEL SECTION---
echo "Sending Parallel Jobs"
for M in "${MATRICI[@]}"; do
    M_NAME=$(basename "$M" .mtx)
    for S in "${SCHEDULERS[@]}"; do
        for T in "${THREADS[@]}"; do
            
            JOB_NAME="${M_NAME}_${S}_${T}"
            OUT_FILE="${RESULTS_DIR}/${JOB_NAME}.out"
            ERR_FILE="${RESULTS_DIR}/${JOB_NAME}.err"
	    
	    check_queue_limit
            
	    echo "Sending: $JOB_NAME"
            
            # Sending Job
            qsub -N "$JOB_NAME" \
                 -o "$OUT_FILE" \
                 -e "$ERR_FILE" \
                 -l select=1:ncpus=$T:mem=$MEMORY \
                 -l walltime=$WALLTIME \
                 -v ROOT="$ROOT",MATRICE="$M",SCHEDULER="$S",NCPUS=$T \
                 "$PBS_TEMPLATE"
	    sleep 1
        done
    done
done

echo "BATCH COMPLETED!"

