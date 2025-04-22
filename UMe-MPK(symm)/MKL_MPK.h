
#ifndef MKL_H
#define MKL_H

#include <mkl_spblas.h>
#include <omp.h>

using namespace std;

void MKL_MPK(int            m,
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

    double start = omp_get_wtime();


    sparse_status_t stat;

    matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.mode = SPARSE_FILL_MODE_LOWER;
    descr.diag = SPARSE_DIAG_NON_UNIT;

    double alpha = 1.0;
    double beta = 0.0;
    
    int *row_ptr_b = new int[m];
    int *row_ptr_e = new int[m];
    for(int i = 0; i < m; i++){
        row_ptr_b[i] = csrRowPtrA[i];
        row_ptr_e[i] = csrRowPtrA[i + 1];
    }

    sparse_matrix_t A;
    stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, m, m, row_ptr_b, row_ptr_e, csrColIdxA, csrValA);

    descr.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descr.mode = SPARSE_FILL_MODE_UPPER;
    descr.diag = SPARSE_DIAG_NON_UNIT;
    
    mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, BENCH_REPEAT+1);
    mkl_sparse_optimize(A);

    if(stat != SPARSE_STATUS_SUCCESS){
        cout << "create failed" << endl;
    }


    double end = omp_get_wtime();
	*pre_time = (end - start) * 1000;



    VALUE_TYPE *x_ref = new VALUE_TYPE[m];
    #pragma omp parallel for
    for(int i = 0; i < m; i++)
        x_ref[i] = x[i];

    // warmup
    stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y);

    *solve_time = 0;
    for(int tt = 0; tt < BENCH_REPEAT; tt++)
    {
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            x[i] = x_ref[i];

        start = omp_get_wtime();

        for(int i = 0; i + 1 < K; i+=2)
        {   
            stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y);
            stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, y, beta, x);
        }
        if (K&1)
            stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A, descr, x, beta, y);
        
        end = omp_get_wtime();
        *solve_time += (end - start) * 1000;
    }

    if (K%2==0)
        #pragma omp parallel for
        for(int j = 0; j < m; j++) y[j] = x[j];

    *solve_time /= BENCH_REPEAT;
    *gflops = 2ll * nnzA * K / (1e6 * (*solve_time));

    delete[] row_ptr_b;
    delete[] row_ptr_e;
    delete[] x_ref;
}

#endif
