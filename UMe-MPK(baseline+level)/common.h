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
#define BENCH_REPEAT 20
#endif



#ifndef CHECK_CORRECTNESS
#define CHECK_CORRECTNESS

bool CheckCorrectness(int m, VALUE_TYPE *y_ref, VALUE_TYPE *y)
{/**************************************CheckCorrectness******************************************/
    
    int error_count = 0;
    for (int i = 0; i < m; i++)
        //both relative error and absolute error should be considered.
        if (fabs(y_ref[i] - y[i]) > 1e-5 * max(fabs(y_ref[i]), fabs(y[i])) && fabs(y_ref[i] - y[i]) > 1e-5)
        {
            error_count++;
            // if (error_count <= 100)
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

void reorder_matrix(int m, int *new_rowPtr, int *new_colIdx, VALUE_TYPE *new_val, int *csrRowPtr,
                                 int *csrColIdx, VALUE_TYPE *csrVal, int *ordering, int *reverse)
{    
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        int row = ordering[i];
        new_rowPtr[i] = j;
        j += csrRowPtr[row + 1] - csrRowPtr[row];
    }
    new_rowPtr[m] = csrRowPtr[m];

#pragma omp parallel for
    for (int i = 0; i < m; i++)
    {
        int row = ordering[i];
        int row_start = csrRowPtr[row];
        int row_end = csrRowPtr[row + 1];
        int tot = row_end - row_start;
        pair<int,VALUE_TYPE> *pairs = new pair<int,VALUE_TYPE>[tot];
        for (int k = row_start; k < row_end; k++)
        {
            pairs[k-row_start].first = reverse[csrColIdx[k]];
            pairs[k-row_start].second = csrVal[k];
        }
        row_start = new_rowPtr[i];
        sort(pairs, pairs+tot);
        for (int k = 0; k < tot; k++)
        {
            new_colIdx[row_start+k] = pairs[k].first;
            new_val[row_start+k] = pairs[k].second;
        }
        delete[] pairs;
    }
}

void reorder_vector(int m, VALUE_TYPE *&new_x, VALUE_TYPE *&x, int *ordering){
    #pragma omp parallel for
    for(int i = 0; i < m; i++)
        new_x[i] = x[ordering[i]];
}


#endif
