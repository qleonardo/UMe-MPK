#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>

#include <sys/time.h>
#include "DC.h"
#include "omp.h"

using namespace std;

#ifndef VALUE_TYPE
#define VALUE_TYPE double
#endif


#ifndef BENCH_REPEAT
#define BENCH_REPEAT 100
#endif



#ifndef CHECK_CORRECTNESS
#define CHECK_CORRECTNESS

bool CheckCorrectness(int m, VALUE_TYPE *y_ref, VALUE_TYPE *y)
{/**************************************CheckCorrectness******************************************/
    
    int error_count = 0;
    for (int i = 0; i < m; i++)
        //both relative error and absolute error should be considered.
        if (fabs(y_ref[i] - y[i]) > 1e-4 * max(fabs(y_ref[i]), fabs(y[i])) && fabs(y_ref[i] - y[i]) > 1e-4)
        {
            error_count++;
            // cout << "ROW [ " << i << " ], NNZ SPAN: "
            //         << "\t ref = " << y[i]
            //         << ", \t std = " << y_ref[i]
            //     << ", \t error = " << fabs(y_ref[i] - y[i])
            //         << endl << flush;
        //    break;
        }

    if (error_count == 0)
    {
        cout << "Check... PASS! \n" << flush;
        return 1;
    }
    else
    {
        cout << "Check... Error! Error_coumnt: " << error_count << "\n" << flush;
        return 0;
    }
}

#endif



#ifndef REORDERING
#define REORDERING

void reorder_matrix(int m, int *new_sssRowPtrA, int *new_sssColIdxA, VALUE_TYPE *new_sssValA, VALUE_TYPE *new_sssDValA, int *sssRowPtrA, int *sssColIdxA, VALUE_TYPE *sssValA, VALUE_TYPE *sssDValA, int *RowRev, int *RowPerm)
{
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        int dst = RowRev[i];
        new_sssRowPtrA[i] = j;
        j += sssRowPtrA[dst + 1]-sssRowPtrA[dst];
    }
    new_sssRowPtrA[m] = sssRowPtrA[m];

#pragma omp parallel for
    for (int i = 0; i < m; i++)
    {
        int dst = RowRev[i];
        int row_start = sssRowPtrA[dst];
        int row_end = sssRowPtrA[dst + 1];
        int tot = row_end -row_start;
        new_sssDValA[i] = sssDValA[dst];
        pair<int,VALUE_TYPE> *pairs = new pair<int,VALUE_TYPE>[tot];
        for (int k = row_start; k < row_end; k++)
        {
            pairs[k-row_start].first = RowPerm[sssColIdxA[k]];
            pairs[k-row_start].second = sssValA[k];
        }
        row_start = new_sssRowPtrA[i];
        sort(pairs, pairs+tot);
        for (int k = 0; k < tot; k++)
        {
            new_sssColIdxA[row_start+k] = pairs[k].first;
            new_sssValA[row_start+k] = pairs[k].second;
        }
        delete[] pairs; 
    }
}

void reorder_vector(int m, VALUE_TYPE *&new_x, VALUE_TYPE *&x, int *ordering){
    #pragma omp parallel for
    for(int i = 0; i < m; i++)
        // new_x[i*(K+1)] = x[ordering[i]];
        new_x[i] = x[ordering[i]];
}


// void matrix_split(int m, int *rowPtrA, int *colIdxA, VALUE_TYPE *valA, int                  *rowPtrL, int *colIdxL, VALUE_TYPE *valL, int *rowPtrU, int *colIdxU, VALUE_TYPE *valU)
// {    
//     int nnzL = 0, nnzU = 0;
//     for(int i = 0; i < m; i++)
//     {
//         bool flag = 0;
//         rowPtrL[i] = nnzL;
//         rowPtrU[i] = nnzU;
//         for(int j = rowPtrA[i]; j < rowPtrA[i+1]; j++)
//         {
//             if (colIdxA[j] < i)
//             {
//                 colIdxL[nnzL] = colIdxA[j];
//                 valL[nnzL++] = valA[j];
//             }
//             else if (colIdxA[j] > i)
//             {
//                 if (!flag)
//                 {
//                     colIdxL[nnzL] = colIdxU[nnzU] = i;
//                     valL[nnzL++] = valU[nnzU++] = 0;
//                     flag = 1;
//                 }
//                 colIdxU[nnzU] = colIdxA[j];
//                 valU[nnzU++] = valA[j];
//             }
//             else
//             {
//                 colIdxL[nnzL] = colIdxU[nnzU] = colIdxA[j];
//                 valL[nnzL++] = valU[nnzU++] = valA[j];
//                 flag = 1;
//             }
//         }
//         if (!flag)
//         {
//             colIdxL[nnzL] = colIdxU[nnzU] = i;
//             valL[nnzL++] = valU[nnzU++] = 0;
//             flag = 1;
//         }
//     }
//     rowPtrL[m] = nnzL;
//     rowPtrU[m] = nnzU;
// }

#endif



#ifndef SSS_TO_CSR
#define SSS_TO_CSR
void convert_SSS_to_CSR(int m, int *sssRowPtrA, int *sssColIdxA, VALUE_TYPE *sssValA, VALUE_TYPE *sssDValA, int *&csrRowPtrA, int *&csrColIdxA, VALUE_TYPE *&csrValA, int start_CSR)
{

    int threads = omp_get_max_threads();
    vector<pair<int,VALUE_TYPE> > *val[threads];
#pragma omp parallel for
    for(int i = 0; i < threads; i++)
        val[i] = new vector<pair<int,VALUE_TYPE> >[m];
#pragma omp parallel for
    for(int i = start_CSR; i < m; i++)
    {
        int id = omp_get_thread_num(); 
        for(int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++)
        {
            val[id][i].push_back(make_pair(sssColIdxA[j], sssValA[j]));
            val[id][sssColIdxA[j]].push_back(make_pair(i, sssValA[j]));
        }
        if (fabs(sssDValA[i])>1e-6)
            val[id][i].push_back(make_pair(i, sssDValA[i]));
    }

    int nnz = 0;
    vector<pair<int, VALUE_TYPE> > *Val = new vector<pair<int, VALUE_TYPE> >[m];
#pragma omp parallel for reduction(+:nnz)
    for (int i = 0; i < m; i++)
    {
        for(int j = 0; j < threads; j++)
        {
            Val[i].insert(Val[i].end(), val[j][i].begin(), val[j][i].end());
            val[j][i].clear();
        }
        sort(Val[i].begin(), Val[i].end());
        nnz += Val[i].size();
    }
#pragma omp parallel for
    for(int j = 0; j < threads; j++)
        delete[] val[j];

    
    int tot = 0;
    csrRowPtrA = new int[m+1];
    csrColIdxA = new int[nnz];
    csrValA = new VALUE_TYPE[nnz];
    for (int i = 0; i < m; i++)
    {
        csrRowPtrA[i] = tot;
        tot += Val[i].size();
    }
    csrRowPtrA[m] = tot;
    
// #pragma omp parallel for
    for (int i = 0; i < m; i++)
    {
        int offset = csrRowPtrA[i];
        for(int j = 0; j < Val[i].size(); j++)
        {
            csrColIdxA[offset+j] = Val[i][j].first;
            csrValA[offset+j] = Val[i][j].second;
            // cout << i << " " << Val[i][j].first << " " << Val[i][j].second << endl;
        }
        Val[i].clear();
    }
    delete[] Val;
}
#endif