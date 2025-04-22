
#ifndef BASE_H
#define BASE_H

#include <omp.h>

using namespace std;

void Base_MPK(int            m,
                 int            nnzA,
                 int            *csrRowPtrA,
                 int            *csrColIdxA,
                 VALUE_TYPE     *csrValA,
                 VALUE_TYPE     *x,
                 VALUE_TYPE     *&y,
                 int            K,
                 double         *solve_time,
                 double         *gflops,
                 double         *pre_time)
{

    double start, end;



    VALUE_TYPE *x_ref = new VALUE_TYPE[m];
    #pragma omp parallel for
    for(int i = 0; i < m; i++)
        x_ref[i] = x[i];

    *solve_time = *pre_time = 0;
    for(int tt = 0; tt < BENCH_REPEAT; tt++)
    {
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            x[i] = x_ref[i];

        start = omp_get_wtime();

        
        for(int j = 0; j+1 < K; j += 2)
        {
            #pragma omp parallel for
            for (int i = 0; i < m; i++)
            {
                VALUE_TYPE sum = 0;
                for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
                    sum += x[csrColIdxA[j]] * csrValA[j];
                y[i] = sum;
            }

            #pragma omp parallel for
            for (int i = 0; i < m; i++)
            {
                VALUE_TYPE sum = 0;
                for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
                    sum += y[csrColIdxA[j]] * csrValA[j];
                x[i] = sum;
            }
        }

        if (K&1)
            #pragma omp parallel for
            for (int i = 0; i < m; i++)
            {
                VALUE_TYPE sum = 0;
                for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
                    sum += x[csrColIdxA[j]] * csrValA[j];
                y[i] = sum;
            }

        end = omp_get_wtime();
        *solve_time += (end - start) * 1000;
    }

    if (K%2==0)
        #pragma omp parallel for
        for (int i = 0; i < m; i++)
            y[i] = x[i];

    *solve_time /= BENCH_REPEAT;
    *gflops = 2ll * nnzA * K / (1e6 * (*solve_time));

    delete[] x_ref;
}

#endif
