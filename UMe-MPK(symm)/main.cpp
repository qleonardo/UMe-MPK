#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "read_mtx.h"
#include "Base_MPK.h"
#include "MKL_MPK.h"
#include "UMe_SymMPK.h"
// #include "UMe_FB_MPK.h"
// #include "UMe_FB_BtoB_MPK.h"
#include "omp.h"
using namespace std;

const int L2_CacheSize = 1024 * 1024;
double coeffs[5] = {0.75, 0.8, 0.85, 1.15, 1.2};//{0.9, 0.95, 1.0, 1.05, 1.1};


int main(int argc, char ** argv)
{
    int m, n, nnzA, isSymmetricA;
    int* csrRowPtrA;
    int* csrColIdxA;
    VALUE_TYPE* csrValA;
    
    int argi = 1;
    int threads = omp_get_max_threads();

    int K;
    char  *filename, *matrix_name;
    if(argc > argi)
    {
        filename = argv[argi++];
        matrix_name = argv[argi++];
        K = atoi(argv[argi]);
    }

    // srand(time(NULL));
    int retcode = read_mtx(filename, &m, &n, &nnzA, &csrRowPtrA, &csrColIdxA, &csrValA);
    if (retcode == -5)
    {
        ofstream performance_result;
        performance_result.open("performance_result.csv", ios::app);
        performance_result << matrix_name << "," << m << "," << nnzA << ",";
        performance_result << endl;
        performance_result.close();
    }
    
    if (retcode < 0) return 0;

    // m=n=5;
    // nnzA=13;
    // int csrRowPtrA[]={0, 2, 4, 6, 11, 13};
    // int csrColIdxA[]={0, 3, 1, 3, 2, 3, 0, 1, 2, 3, 4, 3, 4};
    // VALUE_TYPE csrValA[]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    printf("input matrix %s: ( %i, %i ) nnz = %i\n\n", matrix_name, m, n, nnzA);
    cout << flush;

    
    ofstream performance_result;
    performance_result.open("performance_result.csv", ios::app);
    performance_result << matrix_name << "," << m << "," << nnzA << ",";


    VALUE_TYPE *x = new VALUE_TYPE[m];
    VALUE_TYPE *y = new VALUE_TYPE[m];
    VALUE_TYPE *x_ref = new VALUE_TYPE[m];
    VALUE_TYPE *y_ref = new VALUE_TYPE[m];
    #pragma omp parallel for
    for (int i = 0; i < m; i++) 
        x_ref[i] = x[i] = (rand() % 100 - 100) / 100.0;
        // x_ref[i] = x[i] = 1.0;


    double Base_solve_time, Base_gflops, Base_pre_time;
    Base_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y_ref, K, &Base_solve_time, &Base_gflops, &Base_pre_time);
    printf("Base_solve_time: %.5f\tBase_gflops: %.5f\tBase_pre_time: %.5f\n\n\n", Base_solve_time, Base_gflops, Base_pre_time);
    cout << flush;
    // performance_result << Base_gflops << ",";



    #pragma omp parallel for
    for (int i = 0; i < m; i++) 
        x[i] = x_ref[i];

    double MKL_solve_time, MKL_gflops, MKL_pre_time;
    MKL_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &MKL_solve_time, &MKL_gflops, &MKL_pre_time);
    CheckCorrectness(m, y_ref, y);
    printf("MKL_solve_time: %.5f\tMKL_gflops: %.5f\tMKL_pre_time: %.5f\n\n", MKL_solve_time, MKL_gflops, MKL_pre_time);
    cout << flush;
    performance_result << MKL_pre_time << "," << MKL_gflops << ",";





    int Level_final;
    double coeff_final;
    double UMe_Sym_solve_time_final, UMe_Sym_gflops_final = 0, UMe_Sym_pre_time_final;
    int partSize;
    for(int j = 0; j < 5; j++)
        for(int Level = 2; Level <= 6; Level++)
        {
            partSize = L2_CacheSize * coeffs[j] / 12;
            #pragma omp parallel for
            for (int i = 0; i < m; i++) 
                x[i] = x_ref[i];
                
            double UMe_Sym_solve_time, UMe_Sym_gflops, UMe_Sym_pre_time;
            UMe_Sym_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &UMe_Sym_solve_time, &UMe_Sym_gflops, &UMe_Sym_pre_time, partSize, Level);
            CheckCorrectness(m, y_ref, y);
            printf("UMe_Sym_solve_time: %.5f\tUMe_Sym_gflops: %.5f\tUMe_Sym_pre_time: %.5f\n", UMe_Sym_solve_time, UMe_Sym_gflops, UMe_Sym_pre_time);
            
            if (UMe_Sym_gflops > UMe_Sym_gflops_final)
            {
                coeff_final = coeffs[j];
                Level_final = Level;
                UMe_Sym_solve_time_final = UMe_Sym_solve_time;
                UMe_Sym_gflops_final = UMe_Sym_gflops;
                UMe_Sym_pre_time_final = UMe_Sym_pre_time;
            }
        }



    printf("UMe_Sym_solve_time: %.5f\tUMe_Sym_gflops: %.5f\tUMe_Sym_pre_time: %.5f\n\n", UMe_Sym_solve_time_final, UMe_Sym_gflops_final, UMe_Sym_pre_time_final);
    cout << flush;

    performance_result << coeff_final << ","  << Level_final << ","  << UMe_Sym_pre_time_final << "," << UMe_Sym_gflops_final << ",";




    // performance_result << Level_final << ","  << partSize_final << ","  << UMe_FB_BtoB_pre_time_final << "," << UMe_FB_BtoB_gflops_final;




    performance_result << endl;
    performance_result.close();
    puts("");

    delete[] csrRowPtrA;
    delete[] csrColIdxA;
    delete[] csrValA;

    delete[] x;
    delete[] y;
    delete[] x_ref;
    delete[] y_ref;
}

