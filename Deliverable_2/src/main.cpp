#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include "mmio.h"

using namespace std;

#define MAXSIZE 20.0
#define ITERATIONS 50

class element {
public:
    int row_numb;
    int column_numb;
    double value;

    element() : row_numb(0), column_numb(0), value(0.0) {}
    element(int r, int c, double v) : row_numb(r), column_numb(c), value(v) {}

    bool operator<(const element& x) const {
        if (x.row_numb != row_numb)
        {
            return row_numb < x.row_numb;
        }
        return column_numb < x.column_numb;
    }

    static bool compareForDistribution(const element& a, const element& b, int size) {
        int owner_a = a.row_numb % size;
        int owner_b = b.row_numb % size;

        if (owner_a != owner_b) return owner_a < owner_b;

        return a < b;
    }
};

vector<double> generate_vector(int size){
    vector<double> v(size);
    
    for (int i = 0; i < size; i++)
    {
        v[i] = (double)rand() / MAXSIZE;
    }
    return v;
}


vector<element> readMatrix(const char* filename, int& M, int& N, int& nz){
    MM_typecode matcode;
    FILE *f;

    if ((f = fopen(filename, "r")) == NULL){
            MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0){
        exit(1);
    }

    /* reseve memory for matrices */

    vector<element> matrix;
    matrix.reserve(nz);


    for (int i=0; i<nz; i++)
    {
        int r,c;
        double num;
        fscanf(f, "%d %d %lg\n", &r, &c, &num);
        
        element e(r-1,c-1,num);
        matrix.push_back(e);   
    }

    fclose(f);
	return matrix;
}

vector<element> generateWeakMatrix(int rows_per_proc, int nnz_per_row, int size, int& M, int& N, int& nz) {
    M = rows_per_proc * size;
    N = M;
    nz = M * nnz_per_row;

    vector<element> matrix;
    matrix.reserve(nz);

    // Generazione deterministica basata su seed fisso per riproducibilità
    srand(12); 

    for (int i = 0; i < M; i++) {
        // Elemento diagonale (per evitare righe vuote e garantire stabilità)
        matrix.emplace_back(i, i, 10.0);
        
        // Altri elementi random nella riga
        for (int k = 1; k < nnz_per_row; k++) {
            int col = rand() % N;
            matrix.emplace_back(i, col, ((double)rand() / RAND_MAX));
        }
    }
    return matrix;
}

void spmv_seq(
    int nrow,
    const vector<double>& values,
    const vector<int>& columns,
    const vector<int>& row_ptr,
    const vector<double>& vec_x,
    vector<double>& result)
{
    
    for (int i = 0; i < nrow; i++) {

        double local_sum = 0.0;

        for (int j=row_ptr[i]; j<row_ptr[i+1];j++ ){
            local_sum += values[j] * vec_x[columns[j]];
        }

        result[i]= local_sum;

    }
   
}


int main(int argc, char* argv[]){

    int rank, size;
    int M_global,N_global,nz_global; // M (rows), N (Columns), nz (non-zeros)

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL)+rank);

    vector<element> global_matrix;
    bool weak = false;
    

    if (rank==0){
        if (argc<2){
            if (rank == 0) cerr << "Usage: " << argv[0] << " <filename> OR --weak <rows> <nnz>" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        string arg1 = argv[1];
        if (arg1 == "--weak") {
            if (argc < 4) {
                cerr << "Error: --weak requires rows_per_rank and nnz_per_row" << endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            weak = true;
            int rows_per_rank = atoi(argv[2]);
            int nnz_per_row = atoi(argv[3]);
            global_matrix = generateWeakMatrix(rows_per_rank, nnz_per_row, size, M_global, N_global, nz_global);
        } else {
            global_matrix = readMatrix(argv[1], M_global, N_global, nz_global);
        }

        // Sorting for cyclic distribution
        sort(global_matrix.begin(), global_matrix.end(), [size](const element& a, const element& b) {
            return element::compareForDistribution(a, b, size);
        });

    }


    MPI_Bcast(&M_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz_global, 1, MPI_INT, 0, MPI_COMM_WORLD);


    //-----Data Scatter Preparation------

    //Vectors for the ScatterV
    vector<int> sendcounts(size,0);
    vector<int> displs(size,0);

    //Sorted Vectors for communication
    vector<int> sorted_rows, sorted_cols;
    vector<double> sorted_vals;
    
    if (rank==0){
        
        //We find how many elements we have to send to each process
        for(const auto& e : global_matrix)
        {
            int owner = e.row_numb % size; // 1D Cyclic Distribution
            sendcounts[owner]++;
        }

        //Now we create the displacementfor each process
        int offset=0;
        for (int i = 0; i < size; i++)
        {
            displs[i] = offset;
            offset+=sendcounts[i];
        }
        
        //Convertion from the global_matrix to separate vectors
        sorted_rows.resize(nz_global);
        sorted_cols.resize(nz_global);
        sorted_vals.resize(nz_global);

        for (int i = 0; i < nz_global; i++)
        {
            sorted_rows[i] = global_matrix[i].row_numb / size;
            sorted_cols[i] = global_matrix[i].column_numb;
            sorted_vals[i] = global_matrix[i].value;
        }
        //We free the large Matrix
        vector<element>().swap(global_matrix);
    }

    // --- COMMUNICATION (SCATTERV)---

    int my_nz_count;
    MPI_Scatter(sendcounts.data(), 1, MPI_INT, 
                &my_nz_count, 1, MPI_INT, 
                0, MPI_COMM_WORLD);

    // Receiving Buffers
    vector<int> my_rows(my_nz_count);
    vector<int> my_cols(my_nz_count);
    vector<double> my_vals(my_nz_count);

    // Scatter Vectors
    MPI_Scatterv(sorted_rows.data(), sendcounts.data(), displs.data(), MPI_INT,
                my_rows.data(), my_nz_count, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(sorted_cols.data(), sendcounts.data(), displs.data(), MPI_INT,
                my_cols.data(), my_nz_count, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(sorted_vals.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
                my_vals.data(), my_nz_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    //----CSR CONVERTION----
    
    int my_rows_count = M_global / size;
    if (rank < M_global % size) {
        my_rows_count++;
    }

    //CSR Structure
    vector<double> values = move(my_vals);
    vector<int> columns = move(my_cols);
    vector<int> row_ptr(my_rows_count + 1, 0); 


    for (int r : my_rows) {
        if(r < my_rows_count){
            row_ptr[r + 1]++;
        } 
    }
    // Prefix Sum
    for (int i = 0; i < my_rows_count; i++) {
        row_ptr[i + 1] += row_ptr[i];
    }

    //-----Distributing The X Vector------

    //Now Every process generates a slice of the X vector
    int x_local_size = N_global / size;
    if (rank < N_global % size) {
        x_local_size++;
    }

    vector<double> my_x = generate_vector(x_local_size); 
    vector<double> full_x(N_global);
    vector<double> my_y(my_rows_count, 0.0);

    //Setup AllgatherV
    vector<int> recvcounts_x(size),displs_x(size);
    
    int disp = 0;
    for(int i=0; i<size; i++){
        recvcounts_x[i] = N_global / size;
        if (i < N_global % size) recvcounts_x[i]++;
        displs_x[i] = disp;
        disp += recvcounts_x[i];
    }

    //---- BENCHMARKING ----
    
    // WARMUP
    MPI_Allgatherv(my_x.data(), x_local_size, MPI_DOUBLE, 
                   full_x.data(), recvcounts_x.data(), displs_x.data(), MPI_DOUBLE, MPI_COMM_WORLD);
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
        MPI_Allgatherv(my_x.data(), x_local_size, MPI_DOUBLE, 
                       full_x.data(), recvcounts_x.data(), displs_x.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        t_mid = MPI_Wtime();

        // 2. Computation
        spmv_seq(my_rows_count, values, columns, row_ptr, full_x, my_y);
        t_end = MPI_Wtime();

        local_comm_time += (t_mid - t_start);
        local_comp_time += (t_end - t_mid);
    }

    double t_total_end = MPI_Wtime();
    double local_total = t_total_end - t_total_start;

    // We only consider the maximum time
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

        printf("[CSV_DATA] %s,%d,%d,%e,%e,%e,%.4f,%.2f,%.2f\n", 
               weak ? "WEAK" : "STRONG", 
               M_global,
               size, 
               avg_total, avg_comm, avg_comp, gflops);
        
    }

    MPI_Finalize();
    return 0;
}