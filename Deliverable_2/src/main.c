#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "mmio.h"


#define MAXSIZE 20.0
#define ITERATIONS 50


static int world_size = 1;


typedef struct {
    int row_numb;
    int column_numb;
    double value;
} element;


int compareForDistribution(const void* a, const void* b) {
    const element* element_a = (const element*)a;
    const element* element_b = (const element*)b;

    int owner_a = element_a->row_numb % world_size;
    int owner_b = element_b->row_numb % world_size;

    //Sort on processes
    if (owner_a != owner_b) {
        return owner_a - owner_b;
    }

    // Sort on rows
    if (element_a->row_numb != element_b->row_numb) {
        return element_a->row_numb - element_b->row_numb;
    }

    // Sort on columns
    return element_a->column_numb - element_b->column_numb;
}

// Genera un vettore random
double* generate_vector(int size) {
    double* v = (double*)malloc(size * sizeof(double));
    if (v == NULL) { perror("Malloc Failed!!\n"); exit(1); }

    for (int i = 0; i < size; i++) 
    {
        v[i] = (double)rand() / MAXSIZE;
    }
    return v;
}

// Legge la matrice da file
element* readMatrix(const char* filename, int* M, int* N, int* nz) {
    MM_typecode matcode;
    FILE *f;

    if ((f = fopen(filename, "r")) == NULL) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((mm_read_mtx_crd_size(f, M, N, nz)) != 0) {
        exit(1);
    }

    /* reserve memory for matrices */
    element* matrix = (element*)malloc((*nz) * sizeof(element));
    if (matrix == NULL) { perror("Malloc Failed!!\n"); exit(1); }

    for (int i = 0; i < *nz; i++) 
    {
        int r, c;
        double num;
        fscanf(f, "%d %d %lg\n", &r, &c, &num);
        
        matrix[i].row_numb = r - 1;
        matrix[i].column_numb = c - 1;
        matrix[i].value = num;
    }

    fclose(f);
    return matrix;
}

// Genera matrice weak scaling
element* generateWeakMatrix(int rows_per_proc, int nnz_per_row, int size, int* M, int* N, int* nz) {
    *M = rows_per_proc * size;
    *N = *M;
    *nz = (*M) * nnz_per_row;

    element* matrix = (element*)malloc((*nz) * sizeof(element));
    if (matrix == NULL) { perror("Malloc Failed!!\n"); exit(1); }

    // We use a specific seed for reproducibility
    srand(12); 

    int count = 0;
    for (int i = 0; i < *M; i++) {
        
        //Diagonal Element
        matrix[count].row_numb = i;
        matrix[count].column_numb = i;
        matrix[count].value = 10.0;
        count++;
        
        //Random Elements
        for (int k = 1; k < nnz_per_row; k++) 
        {
            int col = rand() % (*N);
            matrix[count].row_numb = i;
            matrix[count].column_numb = col;
            matrix[count].value = ((double)rand() / RAND_MAX);
            count++;
        }
    }
    return matrix;
}

void spmv_seq(
    int nrow,
    const double* values,
    const int* columns,
    const int* row_ptr,
    const double* vec_x,
    double* result)
{
    for (int i = 0; i < nrow; i++) {
        double local_sum = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) 
        {
            local_sum += values[j] * vec_x[columns[j]];
        }
        result[i] = local_sum;
    }
}

int main(int argc, char* argv[]){

    int rank, size;
    int M_global, N_global, nz_global; // M (rows), N (Columns), nz (non-zeros)

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL)+rank);

    element* global_matrix = NULL;
    int weak = 0; 

    if (rank == 0) {
        if (argc < 2) {
            perror("We are missing some arguments!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (strcmp(argv[1], "--weak") == 0) {
            if (argc < 4) {
                perror("The use of --weak requires rows_per_rank and nnz_per_row\n"); 
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            weak = 1;
            int rows_per_rank = atoi(argv[2]);
            int nnz_per_row = atoi(argv[3]);
            global_matrix = generateWeakMatrix(rows_per_rank, nnz_per_row, size, &M_global, &N_global, &nz_global);
        } else {
            global_matrix = readMatrix(argv[1], &M_global, &N_global, &nz_global);
        }

        // Sorting for cyclic distribution using qsort
        world_size = size;
        qsort(global_matrix, nz_global, sizeof(element), compareForDistribution);
    }

    MPI_Bcast(&M_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz_global, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //-----Data Scatter Preparation------

    //Vectors for the ScatterV
    int* sendcounts = (int*)calloc(size, sizeof(int));
    int* displs = (int*)calloc(size, sizeof(int));

    //Sorted Vectors for communication
    int *sorted_rows = NULL, *sorted_cols = NULL;
    double *sorted_vals = NULL;
    
    if (rank == 0) {
        //We find how many elements we have to send to each process
        for(int i = 0; i < nz_global; i++) 
        {
            int owner = global_matrix[i].row_numb % size; // 1D Cyclic Distribution
            sendcounts[owner]++;
        }

        //Now we create the displacementfor each process
        int offset = 0;
        for (int i = 0; i < size; i++) 
        {
            displs[i] = offset;
            offset+=sendcounts[i];
        }
        
        // Conversione from struct arrays to COO
        sorted_rows = (int*)malloc(nz_global * sizeof(int));
        sorted_cols = (int*)malloc(nz_global * sizeof(int));
        sorted_vals = (double*)malloc(nz_global * sizeof(double));

        for (int i = 0; i < nz_global; i++) 
        {
            sorted_rows[i] = global_matrix[i].row_numb / size;
            sorted_cols[i] = global_matrix[i].column_numb;
            sorted_vals[i] = global_matrix[i].value;
        }

        
        free(global_matrix);
        global_matrix = NULL;
    }

    // --- COMMUNICATION (SCATTERV)---

    int my_nz_count;
    MPI_Scatter(sendcounts, 1, MPI_INT, 
                &my_nz_count, 1, MPI_INT, 
                0, MPI_COMM_WORLD);

    // Receiving Buffers
    int* my_rows = (int*)malloc(my_nz_count * sizeof(int));
    int* my_cols = (int*)malloc(my_nz_count * sizeof(int));
    double* my_vals = (double*)malloc(my_nz_count * sizeof(double));

    // Scatter Vectors
    MPI_Scatterv(sorted_rows, sendcounts, displs, MPI_INT,
                my_rows, my_nz_count, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(sorted_cols, sendcounts, displs, MPI_INT,
                my_cols, my_nz_count, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(sorted_vals, sendcounts, displs, MPI_DOUBLE,
                my_vals, my_nz_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (rank == 0) {
        free(sorted_rows);
        free(sorted_cols);
        free(sorted_vals);
    }
    free(sendcounts);
    free(displs);


    //----CSR CONVERTION----
    
    int my_rows_count = M_global / size;
    if (rank < M_global % size) {
        my_rows_count++;
    }

    // CSR Structure
    double* values = my_vals; 
    int* columns = my_cols;
    int* row_ptr = (int*)calloc(my_rows_count + 1, sizeof(int));

    
    for (int i = 0; i < my_nz_count; i++) 
    {
        int r = my_rows[i];
        if(r < my_rows_count){
            row_ptr[r + 1]++;
        } 
    }

    // Prefix Sum
    for (int i = 0; i < my_rows_count; i++) {
        row_ptr[i + 1] += row_ptr[i];
    }
    
    free(my_rows); 

    //-----Distributing The X Vector------

    //Now Every process generates a slice of the X vector
    int x_local_size = N_global / size;
    if (rank < N_global % size) {
        x_local_size++;
    }

    double* my_x = generate_vector(x_local_size); 
    double* full_x = (double*)malloc(N_global * sizeof(double));
    double* my_y = (double*)calloc(my_rows_count, sizeof(double)); 
    
    // Setup AllgatherV
    int* recvcounts_x = (int*)malloc(size * sizeof(int));
    int* displs_x = (int*)malloc(size * sizeof(int));
    
    int disp = 0;
    for(int i = 0; i < size; i++){
        recvcounts_x[i] = N_global / size;
        if (i < N_global % size) recvcounts_x[i]++;
        displs_x[i] = disp;
        disp += recvcounts_x[i];
    }

    //---- BENCHMARKING ----
    
    // WARMUP
    MPI_Allgatherv(my_x, x_local_size, MPI_DOUBLE, 
                   full_x, recvcounts_x, displs_x, MPI_DOUBLE, MPI_COMM_WORLD);
    spmv_seq(my_rows_count, values, columns, row_ptr, full_x, my_y);

    // Time Variables
    double local_comm_time = 0.0;
    double local_comp_time = 0.0;
    double t_start, t_mid, t_end;

    // Global Synchronization
    MPI_Barrier(MPI_COMM_WORLD); 
    double t_total_start = MPI_Wtime();

    for (int i = 0; i < ITERATIONS; i++) {
        // 1. Comunication
        t_start = MPI_Wtime();
        MPI_Allgatherv(my_x, x_local_size, MPI_DOUBLE, 
                       full_x, recvcounts_x, displs_x, MPI_DOUBLE, MPI_COMM_WORLD);
        t_mid = MPI_Wtime();

        // 2. Computation
        spmv_seq(my_rows_count, values, columns, row_ptr, full_x, my_y);
        t_end = MPI_Wtime();

        local_comm_time += (t_mid - t_start);
        local_comp_time += (t_end - t_mid);
    }

    double t_total_end = MPI_Wtime();
    double local_total = t_total_end - t_total_start;

    // Reduce results
    double max_total, max_comm, max_comp;
    MPI_Reduce(&local_total, &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_comm_time, &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_comp_time, &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double avg_total = max_total / ITERATIONS;
        double avg_comm = max_comm / ITERATIONS;
        double avg_comp = max_comp / ITERATIONS;
        double gflops = (2.0 * nz_global * 1e-9) / avg_total;
        
        printf("GENERAL DEBUG RESULTS:");
        printf("Procs: %d\n", size);
        printf("Matrix: %d x %d (NNZ: %d)\n", M_global, N_global, nz_global);
        printf("Avg Total Time: %e s\n", avg_total);
        printf("  - Comm Time : %e s (%.2f%%)\n", avg_comm, (avg_comm/avg_total)*100);
        printf("  - Comp Time : %e s (%.2f%%)\n", avg_comp, (avg_comp/avg_total)*100);

        printf("[CSV_DATA] %s,%d,%d,%e,%e,%e,%.4f\n", 
               weak ? "WEAK" : "STRONG", 
               M_global,
               size, 
               avg_total, avg_comm, avg_comp, gflops);
    }

    free(values);  
    free(columns); 
    free(row_ptr);
    free(my_x);
    free(full_x);
    free(my_y);
    free(recvcounts_x);
    free(displs_x);

    MPI_Finalize();
    return 0;
}