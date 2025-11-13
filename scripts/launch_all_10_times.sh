#!/bin/bash

# La cartella temporanea che 'launch_all.sh' usa per i risultati
# (deve corrispondere a RESULTS_DIR in launch_all.sh)
CURRENT_RUN_DIR="i_miei_risultati"

# La cartella di destinazione finale dove sposteremo i risultati
# (la stessa che usi nel tuo script Python)
FINAL_RESULTS_BASE_DIR="/home/matteo.berte/results"
mkdir -p "$FINAL_RESULTS_BASE_DIR"

USER="matteo.berte"
USER_A="matteo.b"
echo "--- Avvio 10 iterazioni di benchmark ---"
echo "I risultati finali saranno salvati in $FINAL_RESULTS_BASE_DIR"

# Loop principale: da 1 a 10
for i in $(seq 1 10); do
    echo ""
    echo "================================================="
    echo "=== Avvio Iterazione $i di 10 ==="
    echo "================================================="

    # 1. Pulizia e lancio
    # Pulisce la cartella temporanea per essere sicuri
    rm -rf "$CURRENT_RUN_DIR"
    
    # Lancia la batch di 80 job (o quanti sono)
    echo "Invio batch $i con ./launch_all.sh..."
    ./launch_all.sh
    
    echo "Job inviati. Ora attendo il completamento..."
    echo "(Controllo ogni 60 secondi)"

    # 2. Loop di attesa (Il cuore dell'automazione)
    # Controlla quanti job hai in coda (Q) o in esecuzione (R)
    # 'grep $USER' assicura che contiamo solo i tuoi job
    job_count=$(qstat -u $USER | grep $USER_A | wc -l)

    while [ $job_count -gt 0 ]; do
        echo "$(date '+%H:%M:%S') - In attesa: $job_count job ancora in Q/R..."
        sleep 60
        job_count=$(qstat -u $USER | grep $USER_A | wc -l)
    done

    echo "$(date '+%H:%M:%S') - Batch $i completato! Nessun job in coda."

    # 3. Spostamento e archiviazione
    # Il nome della cartella finale (es. run_1, run_2, ...)
    FINAL_RUN_NAME="run_${i}"
    
    echo "Sposto i risultati da '$CURRENT_RUN_DIR' a '$FINAL_RESULTS_BASE_DIR/$FINAL_RUN_NAME'..."
    mv "$CURRENT_RUN_DIR" "$FINAL_RESULTS_BASE_DIR/$FINAL_RUN_NAME"

    echo "--- Iterazione $i terminata ---"
done

echo ""
echo "================================================="
echo "=== TUTTE LE 10 ITERAZIONI COMPLETATE ==="
echo "================================================="
