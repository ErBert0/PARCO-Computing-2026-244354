#!/bin/bash
#
# --- launch_all.sh ---
# Lancia l'intera suite di benchmark sul cluster PBS.

# --- CONFIGURAZIONE FISSA ---
WALLTIME="03:00:00"
MEMORY=20gb
MAX_JOBS=29
USER="matteo.berte"
USER_A="matteo.b"

echo "--- Avvio batch di benchmark ---"
echo "Memoria richiesta per job: $MEMORY"
echo "Walltime fisso per job:  $WALLTIME"
echo "------------------------------------"


# --- CONFIGURAZIONE BENCHMARK ---
MATRICI_DIR="/home/matteo.berte/Matrices"
MATRICI=(
    "$MATRICI_DIR/1_ASIC.mtx"
    "$MATRICI_DIR/2_kkt.mtx"
    "$MATRICI_DIR/3_mawi.mtx"
    "$MATRICI_DIR/4_nlpkk.mtx"
    "$MATRICI_DIR/5_web-Google.mtx"
)
THREADS=(4 8 16 32 64) # Lista di ncpus da testare
SCHEDULERS=("static" "dynamic" "guided")
RESULTS_DIR="i_miei_risultati" # Cartella per i file .out e .err

# --- ESECUZIONE ---
mkdir -p $RESULTS_DIR


check_queue_limit() {
    # Contiamo quante righe 'matteo.b' appaiono in qstat
    local job_count=$(qstat -u $USER | grep $USER_A |wc -l) 
    
    # Loop di attesa
    while [ $job_count -ge $MAX_JOBS ]; do
        echo "  -> Coda piena ($job_count / $MAX_JOBS). In attesa di 60 secondi..."
        sleep 60
        job_count=$(qstat -u $USER | grep $USER_A | wc -l)
    done
}

# --- 1. BENCHMARK SEQUENZIALE ---
echo "Invio job sequenziali (baseline)..."
for M_PATH in "${MATRICI[@]}"; do
    M_NAME=$(basename "$M_PATH" .mtx)
    
    # Formato nome file: "matrice,scheduler,threads"
    JOB_NAME="${M_NAME}_sequential_1"
    OUT_FILE="${RESULTS_DIR}/${JOB_NAME}.out"
    ERR_FILE="${RESULTS_DIR}/${JOB_NAME}.err"
    
    check_queue_limit

    echo "  - Invio: $JOB_NAME"
    
    qsub -N "$JOB_NAME" \
         -o "$OUT_FILE" \
         -e "$ERR_FILE" \
         -l select=1:ncpus=1:mem=$MEMORY \
         -l walltime=$WALLTIME \
         -v MATRICE="$M_PATH",SCHEDULER="sequential",NCPUS=1 \
         job_template.pbs
    sleep 1
done

# --- 2. BENCHMARK PARALLELO ---
echo "Invio job paralleli..."
for M_PATH in "${MATRICI[@]}"; do
    M_NAME=$(basename "$M_PATH" .mtx)
    for S in "${SCHEDULERS[@]}"; do
        for T in "${THREADS[@]}"; do
            
            JOB_NAME="${M_NAME}_${S}_${T}"
            OUT_FILE="${RESULTS_DIR}/${JOB_NAME}.out"
            ERR_FILE="${RESULTS_DIR}/${JOB_NAME}.err"
	    
	    check_queue_limit
            
	    echo "  - Invio: $JOB_NAME"
            
            # Sottometti il job
            qsub -N "$JOB_NAME" \
                 -o "$OUT_FILE" \
                 -e "$ERR_FILE" \
                 -l select=1:ncpus=$T:mem=$MEMORY \
                 -l walltime=$WALLTIME \
                 -v MATRICE="$M_PATH",SCHEDULER="$S",NCPUS=$T \
                 job_template.pbs
	    sleep 1
        done
    done
done

echo "--- TUTTI I JOB INVIATI ---"

