
#ifndef UMe_LEVEL_MPK_H
#define UMe_LEVEL_MPK_H

#include "DC.h"
#include <omp.h>

using namespace std;

void UMe_LEVEL_MPK(int            m,
                 int            nnzA,
                 int            *csrRowPtrA,
                 int            *csrColIdxA,
                 VALUE_TYPE     *csrValA,
                 VALUE_TYPE     *x,
                 VALUE_TYPE     *&y,
                 int            K,
                 double         *solve_time,
                 double         *gflops,
                 double         *pre_time,
                 int            partSize, 
                 int            Recursion,
                 int            LEVEL)
{
    // Preprocess
    double start = omp_get_wtime();


    //Create DC Tree
    cout << "DCTree begin to create!\n" << flush;
    DC *DCRoot = new DC(m, partSize, Recursion, LEVEL);
    DCRoot->DC_creation(csrColIdxA, csrRowPtrA, m);
    cout << "DCTree created successfully!\n" << flush;

    int *RowRev = DCRoot->DC_get_RowRev();
    int *RowPerm = DCRoot->DC_get_RowPerm();

    // Rerorder
    int *new_rowPtrA = new int[m + 1];
    int *new_colIdxA = new int[nnzA];
    VALUE_TYPE *new_valA = new VALUE_TYPE[nnzA];
    VALUE_TYPE *new_x = new VALUE_TYPE[m];
    reorder_matrix(m, new_rowPtrA, new_colIdxA, new_valA, csrRowPtrA, csrColIdxA, csrValA, RowRev, RowPerm);
    reorder_vector(m, new_x, x, RowRev);

    VALUE_TYPE *y_ref = new VALUE_TYPE[m]();
    char *userArgsA[5] = {(char *)new_rowPtrA, (char *)new_colIdxA, (char *)new_valA, (char *)new_x, (char *)y_ref};
    char *userArgsB[5] = {(char *)new_rowPtrA, (char *)new_colIdxA, (char *)new_valA, (char *)y_ref, (char *)new_x};
    char **userArgs[2] = {userArgsA, userArgsB};


    double end = omp_get_wtime();
	*pre_time = (end - start) * 1000;


    // warm up
    DCRoot->DC_traversal_level(SpMV_kernel_once, SpMV_kernel_twice, userArgs, K);

    *solve_time = 0;
    for(int ii = 0; ii < BENCH_REPEAT; ii++)
    {
        // reset new_x
        reorder_vector(m, new_x, x, RowRev);

        start = omp_get_wtime();
        DCRoot->DC_traversal_level(SpMV_kernel_once, SpMV_kernel_twice, userArgs, K);

        end = omp_get_wtime();
        *solve_time += (end - start) * 1000;
    }

    *solve_time /= BENCH_REPEAT;
    *gflops = 2ll * nnzA * K / (1e6 * (*solve_time));


    if (K%2==0)
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            y[i] = new_x[RowPerm[i]];
    else
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            y[i] = y_ref[RowPerm[i]];



    delete DCRoot;
    delete[] y_ref;
    delete[] new_rowPtrA;
    delete[] new_colIdxA;
    delete[] new_valA;
    delete[] new_x;
}

#endif
