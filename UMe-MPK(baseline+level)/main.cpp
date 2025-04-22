#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "read_mtx.h"
#include "MKL_MPK.h"
#include "Base_MPK.h"
#include "DC.h"
#include "UMe_MPK.h"
#include "UMe_LEVEL_MPK.h"
#include "omp.h"
using namespace std;

const int L2_Cache_Size = 1024 * 1024;
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
    if (read_mtx(filename, &m, &n, &nnzA, &csrRowPtrA, &csrColIdxA, &csrValA) < 0)
	return 0;

    // m=n=5;
    // nnzA=13;
    // int csrRowPtrA[]={0, 2, 4, 6, 11, 13};
    // int csrColIdxA[]={0, 3, 1, 3, 2, 3, 0, 1, 2, 3, 4, 3, 4};
    // VALUE_TYPE csrValA[]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

    printf("input matrix %s: ( %i, %i ) nnz = %i\n", matrix_name, m, n, nnzA);
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


    //std result
    for(int j = 0; j < K; j++)
    {
        #pragma omp parallel for
        for (int i = 0; i < m; i++)
        {
            VALUE_TYPE sum = 0;
            for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
                sum += x[csrColIdxA[j]] * csrValA[j];
            y_ref[i] = sum;
        }
        #pragma omp parallel for
        for (int i = 0; i < m; i++)
            x[i] = y_ref[i];
    }





    #pragma omp parallel for
    for (int i = 0; i < m; i++) 
        x[i] = x_ref[i];

    double Base_solve_time, Base_gflops, Base_pre_time;
    Base_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &Base_solve_time, &Base_gflops, &Base_pre_time);
    CheckCorrectness(m, y_ref, y);
    printf("Base_solve_time: %.5f\tBase_gflops: %.5f\tBase_pre_time: %.5f\n\n", Base_solve_time, Base_gflops, Base_pre_time);
    cout << flush;
    performance_result << Base_pre_time << "," << Base_gflops << ",";





    #pragma omp parallel for
    for (int i = 0; i < m; i++) 
        x[i] = x_ref[i];

    double MKL_solve_time, MKL_gflops, MKL_pre_time;
    MKL_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &MKL_solve_time, &MKL_gflops, &MKL_pre_time);
    CheckCorrectness(m, y_ref, y);
    printf("MKL_solve_time: %.5f\tMKL_gflops: %.5f\tMKL_pre_time: %.5f\n\n", MKL_solve_time, MKL_gflops, MKL_pre_time);
    cout << flush;
    performance_result << MKL_pre_time << "," << MKL_gflops << ",";




    double coeff_final;
    double UMe_solve_time_final, UMe_gflops_final = 0, UMe_pre_time_final;
    for(int j = 0; j < 5; j++)
    {
        double coeff = coeffs[j];
        int partSize = L2_Cache_Size  * coeff / 12;

        #pragma omp parallel for
        for (int i = 0; i < m; i++) 
            x[i] = x_ref[i];
            
        double UMe_solve_time, UMe_gflops, UMe_pre_time;
        UMe_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &UMe_solve_time, &UMe_gflops, &UMe_pre_time, partSize, 2);
        CheckCorrectness(m, y_ref, y);
        printf("UMe_solve_time: %.5f\tUMe_gflops: %.5f\tUMe_pre_time: %.5f\n", UMe_solve_time, UMe_gflops, UMe_pre_time);
        
        if (UMe_gflops > UMe_gflops_final)
        {
            coeff_final = coeff;
            UMe_solve_time_final = UMe_solve_time;
            UMe_gflops_final = UMe_gflops;
            UMe_pre_time_final = UMe_pre_time;
        }
    }

    printf("UMe_solve_time: %.5f\tUMe_gflops: %.5f\tUMe_pre_time: %.5f\n\n\n\n", UMe_solve_time_final, UMe_gflops_final, UMe_pre_time_final);
    cout << flush;
    performance_result << coeff_final << ","  << UMe_pre_time_final << "," << UMe_gflops_final << ",";
	
    


    int Level_final;
    double UMe_LEVEL_solve_time_final, UMe_LEVEL_gflops_final = 0, UMe_LEVEL_pre_time_final;
    for(int j = 0; j < 5; j++)
        for(int Level = max(5, K/2); Level <= K; Level++)
        {
            double coeff = coeffs[j];
            int partSize = L2_Cache_Size  * coeff / 12;

            #pragma omp parallel for
            for (int i = 0; i < m; i++) 
                x[i] = x_ref[i];
                
            double UMe_LEVEL_solve_time, UMe_LEVEL_gflops, UMe_LEVEL_pre_time;
            UMe_LEVEL_MPK(m, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, K, &UMe_LEVEL_solve_time, &UMe_LEVEL_gflops, &UMe_LEVEL_pre_time, partSize, 2, Level);
            CheckCorrectness(m, y_ref, y);
            printf("UMe_LEVEL_solve_time: %.5f\tUMe_LEVEL_gflops: %.5f\tUMe_LEVEL_pre_time: %.5f\n", UMe_LEVEL_solve_time, UMe_LEVEL_gflops, UMe_LEVEL_pre_time);
            
            if (UMe_LEVEL_gflops > UMe_LEVEL_gflops_final)
            {
                Level_final = Level;
                coeff_final = coeff;
                UMe_LEVEL_solve_time_final = UMe_LEVEL_solve_time;
                UMe_LEVEL_gflops_final = UMe_LEVEL_gflops;
                UMe_LEVEL_pre_time_final = UMe_LEVEL_pre_time;
            }
        }

    printf("UMe_LEVEL_solve_time: %.5f\tUMe_LEVEL_gflops: %.5f\tUMe_LEVEL_pre_time: %.5f\n\n\n\n", UMe_LEVEL_solve_time_final, UMe_LEVEL_gflops_final, UMe_LEVEL_pre_time_final);
    cout << flush;

    performance_result << coeff_final << "," << Level_final << ","  << UMe_LEVEL_pre_time_final << "," << UMe_LEVEL_gflops_final << ",";




    performance_result << "\n";
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

