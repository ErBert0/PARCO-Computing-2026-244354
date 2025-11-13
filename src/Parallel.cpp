#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <sstream>
#include <omp.h>
using namespace std;

#define MAXSIZE 20

class element{

public:
    double value;
    int row_numb;
    int column_numb;

    element(int row, int column,double value){
        this->row_numb = row;
        this->column_numb = column;
        this->value = value;
    }

    //Per ordinare la matrice
    bool operator<(const element& x) const {

        if (x.row_numb != row_numb)
        {
            return row_numb < x.row_numb;
        }
        return column_numb < x.column_numb;
    }

private:

    void set_value(int num){
        this->value = num;
    }

    void set_row(int num){
        this->row_numb = num;
    }

    void set_column(int num){
        this->column_numb = num;
    }

};


vector<double> generate_vector(int size){
    vector<double> vector(size);
    
    for (int i = 0; i < size; i++)
    {
        vector[i] = (double)rand() / MAXSIZE;
    }
    return vector;
}

vector<double> spmv_parallel(
    const vector<double>& vec,
    int nrow,
    const vector<double>& values,
    const vector<int>& columns,
    const vector<int>& row_ptr,
    string scheduler)
{
    
    vector<double> result(nrow);


    if (scheduler == "static")
    {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nrow; i++) {

        double local_sum = 0;

        for (int j=row_ptr[i]; j<row_ptr[i+1];j++ ){
            local_sum += values[j] * vec[columns[j]];
        }

        result[i]= local_sum;
        }
    }
    else if (scheduler == "dynamic")
    {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nrow; i++) {

        double local_sum = 0;

        for (int j=row_ptr[i]; j<row_ptr[i+1];j++ ){
            local_sum += values[j] * vec[columns[j]];
        }

        result[i]= local_sum;
        }
    }
    else if (scheduler == "guided")
    {
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < nrow; i++) {

        double local_sum = 0;

        for (int j=row_ptr[i]; j<row_ptr[i+1];j++ ){
            local_sum += values[j] * vec[columns[j]];
        }

        result[i]= local_sum;
        }
    }

    return result;
        
}


int main(int argc, char* argv[]){

    srand(time(NULL));
    

    if (argc != 3)
    {
        cout << "Ci sono dei parametri mancanti"<<endl;
        exit(1);
    }

    string matrix_name = argv[1];
    string scheduler = argv[2];

    string row;
    int n_row, n_col, n_nonnull;


    //APERTURA E LETTURA FILE
    ifstream matrixfile(matrix_name);

    if(!matrixfile.is_open()){
        cout << "Errore apertura File" << endl;
    }

    while(getline(matrixfile, row )){
        
        if (row[0] == '%'){
            continue;
        }
        
        istringstream iss(row); //Crea uno stream della stringa come se fosse un file
        iss >> n_row >> n_col >> n_nonnull;
        break;
    }

    //---Creazione matrice---
    vector<element> matrice;
    
    while (getline(matrixfile, row))
    {
        istringstream iss(row);
        int row,col;
        double num;
        iss >> row >> col>> num;

        element e(row-1,col-1,num);
        matrice.push_back(e);
    }

    // ---Ordinamento Matrice---
    sort(matrice.begin(), matrice.end());


    //Conversione in CSR
    vector<double> values;
    vector<int> columns;
    vector<int> row_ptr(n_row+1);

    row_ptr[0] = 0;
    int current_row = 0;

    for (const auto &e : matrice)
    {
        values.push_back(e.value);
        columns.push_back(e.column_numb);

        while (current_row < e.row_numb)
        {
            row_ptr[current_row+1] = values.size() -1;
            current_row++; 
        }
        
    }


    while (current_row < n_row)
    {
        row_ptr[current_row+1] = n_nonnull;
        current_row++;
    }


    // --- 1. Genera Vettore Random ---
    vector<double> x = generate_vector(n_col); 


    // --- 2. Inizio Testing ---
    struct timespec t0m, t1m;

    clock_t t0 = clock();
    clock_gettime(CLOCK_MONOTONIC, &t0m);
    vector<double> y_par = spmv_parallel(x, n_row, values, columns, row_ptr, scheduler);
    clock_t t1 = clock();
    clock_gettime(CLOCK_MONOTONIC, &t1m);

    double elapsed_ms = (t1m.tv_sec - t0m.tv_sec)*1000.0 +
                     (t1m.tv_nsec - t0m.tv_nsec) / 1e6;

    
    // Stampa il tempo CPU in secondi (va bene così)
    cout << "CPU time:" << static_cast<double>(t1 - t0) / CLOCKS_PER_SEC << " s";

    // Stampa il wall-time in millisecondi
    cout << " | Elapsed time:" << elapsed_ms << " ms" << endl;


    return 0;
}