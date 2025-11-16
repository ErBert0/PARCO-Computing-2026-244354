# Parallel Computing Project: SpMV Algorithm


## 1. Structure

- **/src**: Contains the C++ source code.
- **/scripts**: Contains all utility scripts for running benchmarks on the cluster.
- **/results**: Contains the final aggregated output files (`.xlsx`) with all the data collected.
- **/plots**: Contains the final graphs and figures used in the report.

## 2. Compiler Version and Environment

All benchmarks were conducted on the university cluster, using:

- **Compiler:** `gcc91` (GCC 9.1)
- **Compilation Flags:** `g++-9.1.0 -O3 -fopenmp -o <output_name> <source_file>`

## 3. How to Compile and Run

### Step 0: Clone the Repository
```bash
git clone https://github.com/ErBert0/PARCO-Computing-2026-244354.git
cd PARCO-Computing-2026-244354
```

### Step 1: Compilation

To compile the `Sequential.cpp` and `Parallel.cpp` source codes, you can run the provided helper script:

```bash
# This script will compile both .cpp files from /src
# Here you will be albe to change the compiler flags or the source paths.
./scripts/compile.sh
``` 
This will create `Sequential.out` and `Parallel.out` in the main directory.

### Step 2: Configure Paths

You will have to manually edit **two different variables** in the `script/launch_all.sh` file:
- `MATRICI_DIR`: Must Point to the absolute path of your matrices folder on the cluster
- `MATRICI`: Must contain the correct names of all your matrices

### Step 3: Running the Benchmark

To run the entire benchmark you can use the main script called `launch_all_10_times.sh`
- **Important:** I suggest using a terminal multiplexer such as "screen". The script may take several hours to complete.
- To launch the main script simply run:
    ```bash
    ./scripts/launch_all_10_times.sh
    ```
- The script will create 10 new directories in `/results` named `run_1`,`run_2` etc... Each containing the `.out` and `.err` files for that specific run.

### Step 4 (Optional): Modifying parameters & manual Runs

If you wish to modify the benchmark parameters you will have to edit the `script/launch_all.sh` file. Here in the top section you will find all the main variables.

- **Matrices:** To change the list of matrices, edit the `MATRICI` array.
- **Threads:** To change the thread counts, edit the `THREADS` array .
- **Memory/Walltime:** These job parameters can be modified at the top of the `scripts/launch_all.sh` file.

To **submit a job manually:**
```bash
qsub -N "manual_submit" \
     -o "test.out" \
     -e "test.err" \
     -l select=1:ncpus=16:mem=20gb \
     -l walltime=03:00:00 \
     -v MATRICE="<path_to_the_matrix>/1_ASIC.mtx",SCHEDULER="dynamic",NCPUS=16 \
     scripts/job_template.pbs
```

## 4. Input and Output

- **Input:** The C++ program accepts the path to a matrix in `.mtx` format. The parallel version also accepts a scheduler name (`static`, `dynamic`, or `guided`) as the second argument.
- **Output:** You will end up with two different files `.out` and `.err`
    - The `.out` file will contain the CPU time and the elapsed time of the process in ms
    - The `.err` file will contain the `perf stat` analysis 


