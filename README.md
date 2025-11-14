# Progetto di Calcolo Parallelo: Ottimizzazione SpMV

[cite_start]**Autore:** Mario Rossi (Sostituisci con nome, matricola, email) [cite: 128]
[cite_start]**Link alla Repository:** (Link alla tua repository Git) [cite: 128]
[cite_start]**Corso:** Introduction to Parallel Computing (a.y. 2025/2026) [cite: 1, 6]

---

## 1. Descrizione del Progetto

Questo progetto implementa e valuta le performance di una moltiplicazione matrice-vettore per matrici sparse (SpMV) in C++ utilizzando OpenMP.

L'obiettivo è analizzare come diverse strategie di scheduling (`static`, `dynamic`, `guided`) e un numero variabile di thread influenzino le performance, con un'analisi approfondita dei bottleneck di memoria tramite `perf`.

## 2. Struttura della Repository

[cite_start]La repository è strutturata come segue, seguendo le linee guida del corso [cite: 145-156]:

- [cite_start]**/src**: Contiene il codice sorgente C++ (`Sequential.cpp`, `Parallel.cpp`)[cite: 150].
- [cite_start]**/scripts**: Contiene tutti gli script per l'esecuzione dei benchmark sul cluster e l'analisi dei risultati[cite: 152].
- [cite_start]**/results**: Contiene i file di output aggregati (`.xlsx`) generati dagli script di parsing[cite: 154].
- [cite_start]**/plots**: Contiene i grafici e le figure finali utilizzate nella relazione[cite: 156].

## 3. Ambiente e Toolchain

Le analisi sono state condotte sull'ambiente cluster (specifica il nome del cluster).

- [cite_start]**Compilatore:** `gcc91` (GCC 9.1) [cite: 158, 162]
- [cite_start]**Flag di Compilazione:** `g++ -O3 -fopenmp -o <nome_output> <file_sorgente>` [cite: 158]

## 4. Come Compilare ed Eseguire

[cite_start]Seguire questi passaggi per riprodurre gli esperimenti sul cluster[cite: 159, 211].

### A. Compilazione

1.  Caricare il modulo GCC corretto:
    ```bash
    module load gcc91
    ```
2.  Compilare il codice sequenziale e parallelo nella cartella `src/`:
    ```bash
    g++ -O3 -fopenmp -o Sequential.out src/Sequential.cpp
    g++ -O3 -fopenmp -o Parallel.out src/Parallel.cpp
    ```

### B. Esecuzione dei Benchmark (Completa)

Per eseguire l'intera suite di test (10 run per 80+ configurazioni):

1.  Assicurati che gli script `launch_all.sh` e `run_all_10_iterations.sh` abbiano i permessi di esecuzione:
    ```bash
    chmod +x scripts/launch_all.sh
    chmod +x scripts/run_all_10_iterations.sh
    ```
2.  Avvia lo script "master" (consigliato l'uso di `screen` o `tmux`):
    ```bash
    screen -S benchmark_spmv
    ./scripts/run_all_10_iterations.sh
    ```
3.  Lo script creerà automaticamente 10 cartelle (`run_1`, `run_2`, ...) nella cartella `/home/mat/Work/results` (come specificato nello script) e attenderà il completamento di ogni batch.

### C. Esecuzione di un Singolo Test (Manuale)

Per lanciare un singolo job manualmente (es. matrice `1_ASIC` con scheduler `dynamic` e 8 thread):

```bash
qsub -N "test_manuale" \
     -o "test.out" \
     -e "test.err" \
     -l select=1:ncpus=8:mem=20gb \
     -l walltime=03:00:00 \
     -v MATRICE="/home/matteo.berte/Matrices/1_ASIC.mtx",SCHEDULER="dynamic",NCPUS=8 \
     scripts/job_template.pbs
```

## 5. Analisi dei Risultati

Per parsare tutti i file `.out` e `.err` e generare il file Excel finale con i P90 e le statistiche `perf`:

1.  Assicurati che il percorso `general_dir` nello script `Parsing.py` sia corretto.
2.  Esegui lo script Python (richiede `pandas`, `numpy`, `openpyxl`):
    ```bash
    python scripts/Parsing.py
    ```
3.  Questo creerà il file `benchmark_completo_perf.xlsx` nella cartella `scripts/`.

## [cite_start]6. Input, Output e Parametri [cite: 160, 161]

- **Input:** Il programma accetta il percorso di una matrice in formato `.mtx` e (per il parallelo) il nome dello scheduler (`static`, `dynamic`, `guided`).
- **Output:** Il programma stampa sull'output standard (file `.out`) l'Elapsed Time P90 (dopo la correzione) e sull'errore standard (file `.err`) l'output di `perf stat`.
- [cite_start]**Modifica Parametri[cite: 215]:**
    - **Matrici:** Per cambiare le matrici, modificare la lista `MATRICI` nello script `scripts/launch_all.sh`.
    - **Thread:** Per cambiare i thread, modificare la lista `THREADS` nello script `scripts/launch_all.sh`.
    - **Memoria/Walltime:** Questi parametri possono essere modificati all'inizio di `scripts/launch_all.sh`.