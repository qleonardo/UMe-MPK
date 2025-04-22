
#ifndef UMe_MPK_H
#define UMe_MPK_H

#include "DC.h"
#include <omp.h>

using namespace std;

inline void SpMV_kernel_once(char ** userArgs, DCArgs *treeArgs)
{
    int *csrRowPtrA         = (int *)userArgs[0];
    int *csrColIdxA         = (int *)userArgs[1];
    VALUE_TYPE *csrValA     = (VALUE_TYPE *)userArgs[2];
    VALUE_TYPE *x           = (VALUE_TYPE *)userArgs[3];
    VALUE_TYPE *y           = (VALUE_TYPE *)userArgs[4];
   
    int ns = treeArgs->firstRow;
    int ne = treeArgs->lastRow;

    // if (ns >= ne) return;
    
    for (int i = ns; i < ne; i++)
    {
        VALUE_TYPE sum = 0;

        for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
            sum += x[csrColIdxA[j]] * csrValA[j];
        y[i] = sum;
    }
}

inline void SpMV_kernel_twice(char ** userArgs, DCArgs *treeArgs)
{
    int *csrRowPtrA         = (int *)userArgs[0];
    int *csrColIdxA         = (int *)userArgs[1];
    VALUE_TYPE *csrValA     = (VALUE_TYPE *)userArgs[2];
    VALUE_TYPE *x           = (VALUE_TYPE *)userArgs[3];
    VALUE_TYPE *y           = (VALUE_TYPE *)userArgs[4];
   
    int ns = treeArgs->firstRow;
    int ne = treeArgs->lastRow;

    // if (ns >= ne) return;

    for (int i = ns; i < ne; i++)
    {
        VALUE_TYPE sum = 0;
        for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
            sum += x[csrColIdxA[j]] * csrValA[j];
        y[i] = sum;
    }

    for (int i = ns; i < ne; i++)
    {
        VALUE_TYPE sum = 0;
        for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
            sum += y[csrColIdxA[j]] * csrValA[j];
        x[i] = sum;
    }
}


void UMe_MPK(int            m,
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
                 int            Recursion)
{
    // Preprocess
    double start = omp_get_wtime();


    //Create DC Tree
    cout << "DCTree begin to create!\n" << flush;
    DC *DCRoot = new DC(m, partSize, Recursion);
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


    DCRoot->DC_traversal(SpMV_kernel_once, SpMV_kernel_twice, userArgs, K);

    *solve_time = 0;
    for(int ii = 0; ii < BENCH_REPEAT; ii++)
    {
        // reset new_x
        reorder_vector(m, new_x, x, RowRev);

        start = omp_get_wtime();
        DCRoot->DC_traversal(SpMV_kernel_once, SpMV_kernel_twice, userArgs, K);

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
