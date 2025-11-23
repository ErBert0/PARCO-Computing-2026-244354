#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace std;

#define MAXSIZE 20
#define NUM_RUN 50

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

    //For Sorting
    bool operator<(const element& x) const {

        if (x.row_numb != row_numb)
        {
            return row_numb < x.row_numb;
        }
        return column_numb < x.column_numb;
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

void spmv_seq(
    const vector<double>& vec,
    int nrow,
    const vector<double>& values,
    const vector<int>& columns,
    const vector<int>& row_ptr,
    vector<double>& result)
{
    
    
    for (int i = 0; i < nrow; i++) {

        double local_sum = 0.0;

        for (int j=row_ptr[i]; j<row_ptr[i+1];j++ ){
            local_sum += values[j] * vec[columns[j]];
        }

        result[i]= local_sum;

    }

        
}


int main(int argc, char* argv[]){


    srand(time(NULL));
    
    if (argc != 2)
    {
        cout << "The Matrix file must be specified" <<endl;
        exit(1);
    }

    string matrix_name = argv[1];
    
    
    string row;
    int n_row, n_col, n_nonnull;


    ifstream matrixfile(matrix_name);

    if(!matrixfile.is_open()){
        cout << "Errore apertura File" << endl;
        exit(1);
    }

    while(getline(matrixfile, row )){
        
        if (row[0] == '%'){
            continue;
        }
        
        istringstream iss(row); 
        iss >> n_row >> n_col >> n_nonnull;
        break;
    }

    //Matrix Creation
    vector<element> matrice;
    matrice.reserve(n_nonnull);
    
    while (getline(matrixfile, row))
    {
        istringstream iss(row);
        int r,c;
        double num;
        iss >> r >> c >> num;

        element e(r-1,c-1,num);
        matrice.push_back(e);
    }

    //Sorting
    sort(matrice.begin(), matrice.end());

    //CSR Conversion
    vector<double> values;
    vector<int> columns;
    vector<int> row_ptr(n_row+1);

    values.reserve(n_nonnull);
    columns.reserve(n_nonnull);

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
        row_ptr[current_row+1] = values.size();
        current_row++;
    }


    //1. Create Vectors
    vector<double> x = generate_vector(n_col); 
    vector<double> y(n_row,0.0);
    vector<double> timings;
    timings.reserve(NUM_RUN);

    //CACHE WARMUP
    spmv_seq(x, n_row, values, columns, row_ptr,y);

    struct timespec t0m, t1m;

     //2. Start Testing
    for (int i = 0; i< NUM_RUN; i++){
        
        clock_gettime(CLOCK_MONOTONIC, &t0m);
        spmv_seq(x, n_row, values, columns, row_ptr,y);
        clock_gettime(CLOCK_MONOTONIC, &t1m);
    
        
        double elapsed = (t1m.tv_sec - t0m.tv_sec)*1000.0 +
                     (t1m.tv_nsec - t0m.tv_nsec) / 1e6;
        timings.push_back(elapsed);
    
    }

    sort(timings.begin(),timings.end());
    double p90_time = timings[(int)(NUM_RUN * 0.9) - 1];

    cout << "P90_TIME_MS:" << p90_time << endl;
    return 0;
}