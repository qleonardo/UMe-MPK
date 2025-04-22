
#ifndef UMe_SYM_MPK_H
#define UMe_SYM_MPK_H

#include "DC.h"
#include <omp.h>

using namespace std;


inline void SpMV_kernel_once(char ** userArgs, DCArgs *treeArgs)
{
    int *sssRowPtrA         = (int *)userArgs[0];
    int *sssColIdxA         = (int *)userArgs[1];
    VALUE_TYPE *sssValA     = (VALUE_TYPE *)userArgs[2];
    VALUE_TYPE *sssDValA    = (VALUE_TYPE *)userArgs[3];
    VALUE_TYPE *x           = (VALUE_TYPE *)userArgs[4];
    VALUE_TYPE *y           = (VALUE_TYPE *)userArgs[5];
    VALUE_TYPE *z           = (VALUE_TYPE *)userArgs[6];
    int *csrRowPtrA         = (int *)userArgs[7];
    int *csrColIdxA         = (int *)userArgs[8];
    VALUE_TYPE *csrValA     = (VALUE_TYPE *)userArgs[9];
   
    int ns = treeArgs->firstRow;
    int ne = treeArgs->lastRow;

    for (int i = ns; i < ne; i++)
    {
        z[i] = 0;
        VALUE_TYPE sum = x[i] * sssDValA[i];
        for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
            sum += x[csrColIdxA[j]] * csrValA[j];
        for (int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++) 
        {
            sum += x[sssColIdxA[j]] * sssValA[j];
            y[sssColIdxA[j]] += x[i] * sssValA[j];
        }
        y[i] += sum;
    }
}


inline void SpMV_kernel_twice(char ** userArgs, DCArgs *treeArgs)
{
    int *sssRowPtrA         = (int *)userArgs[0];
    int *sssColIdxA         = (int *)userArgs[1];
    VALUE_TYPE *sssValA     = (VALUE_TYPE *)userArgs[2];
    VALUE_TYPE *sssDValA    = (VALUE_TYPE *)userArgs[3];
    VALUE_TYPE *x           = (VALUE_TYPE *)userArgs[4];
    VALUE_TYPE *y           = (VALUE_TYPE *)userArgs[5];
    VALUE_TYPE *z           = (VALUE_TYPE *)userArgs[6];
   
    int ns = treeArgs->firstRow;
    int ne = treeArgs->lastRow;
    
    for (int i = ns; i < ne; i++)
    {
        z[i] = 0;
        VALUE_TYPE sum = x[i] * sssDValA[i];
        for (int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++) 
        {
            sum += x[sssColIdxA[j]] * sssValA[j];
            y[sssColIdxA[j]] += x[i] * sssValA[j];
        }
        y[i] += sum;
    }

    for (int i = ns; i < ne; i++)
    {
        x[i] = 0;
        VALUE_TYPE sum = y[i] * sssDValA[i];
        for (int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++) 
        {
            sum += y[sssColIdxA[j]] * sssValA[j];
            z[sssColIdxA[j]] += y[i] * sssValA[j];
        }
        z[i] += sum;
    }

}

inline void Final_SpMV(char ** userArgs, DCArgs *treeArgs)
{
    VALUE_TYPE *x           = (VALUE_TYPE *)userArgs[4];
    VALUE_TYPE *y           = (VALUE_TYPE *)userArgs[5];
    VALUE_TYPE *z           = (VALUE_TYPE *)userArgs[6];
    int *csrRowPtrA         = (int *)userArgs[7];
    int *csrColIdxA         = (int *)userArgs[8];
    VALUE_TYPE *csrValA     = (VALUE_TYPE *)userArgs[9];
   
    int ns = treeArgs->firstRow;
    int ne = treeArgs->lastRow;

    for (int i = ns; i < ne; i++)
    {
	z[i] = 0;
        VALUE_TYPE sum = 0;
        for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++) 
            sum += x[csrColIdxA[j]] * csrValA[j];
        y[i] += sum;
    }
}


long long create_Conflict_Graph(int* sssRowPtrA, int* sssColIdxA, int *nRowPerRow, int **Row2Row, int m)
{
    double t1 = omp_get_wtime();

    // cout << "create conflict graph 0" << endl << flush;
    int threads = omp_get_max_threads();
    vector<int> *output[threads];
    int **vis = new int*[threads];
    
#pragma omp parallel for
    for(int i = 0; i < threads; i++)
    {
        vis[i] = new int[m];
        output[i] = new vector<int>[m];
        for(int j = 0; j < m; j++)
            vis[i][j] = -1;
    }

    // cout << "create conflict graph 1" << endl << flush;
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < m; i++)
    {
        int id = omp_get_thread_num();
        for(int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++)
            output[id][sssColIdxA[j]].push_back(i);
    }


    // cout << "create conflict graph 2" << endl << flush;
    long long conflicts = 0;
#pragma omp parallel for reduction(+:conflicts)
    for(int i = 0; i < m; i++)
    {
        int sum = 0;
        for(int t = 0; t < threads; t++)
            sum += output[t][i].size();
        conflicts += 1ll*(sum-1)*sum/2;
    }
    if (conflicts > 1e10) return conflicts;
    conflicts = 0;
            
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < m; i++)
    {
        vector<int> confRow;
        int id = omp_get_thread_num();

        for(int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++)
        {
            confRow.push_back(sssColIdxA[j]);
            vis[id][sssColIdxA[j]] = i;
        }

        for(int t = 0; t < threads; t++)
            for (int j = 0; j < output[t][i].size(); j++)
            {
                int row = output[t][i][j];
                if (vis[id][row] != i)
                {
                    confRow.push_back(row);
                    vis[id][row] = i;
                }
            }

        for(int t = 0; t < threads; t++) 
            for(int j = sssRowPtrA[i]; j < sssRowPtrA[i+1]; j++)
            {
                int jj = sssColIdxA[j];
                for (int k = 0; k < output[t][jj].size(); k++)
                {
                    int row = output[t][jj][k];
                    if (vis[id][row] != i)
                    {
                        confRow.push_back(row);
                        vis[id][row] = i;
                    }
                }

            }

        int k = 0;
        Row2Row[i] = new int[confRow.size()];
        for(int j = 0; j < confRow.size(); j++)
            if (confRow[j] != i)
                Row2Row[i][k++] = confRow[j];
        nRowPerRow[i] = k;
        // cout << k << endl << flush;

        confRow.clear();
    }

    for(int i = 0; i < m; i++)
        conflicts += nRowPerRow[i];
    // cout << "create conflict graph 4" << endl << flush;
#pragma omp parallel for
    for(int i = 0; i < threads; i++)
    {
        for(int j = 0; j < m; j++)
            output[i][j].clear();
        delete[] output[i];
        delete[] vis[i];
    }
    delete[] vis;
    
    double t2 = omp_get_wtime();
    cout << "create conflict graph time: " << t2 - t1 << endl << flush;


    return conflicts;
}


void UMe_Sym_MPK(int            m,
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
                 int            numParts)
{
    // Preprocess
    double start = omp_get_wtime();


    int num = 0;
    int *RowValue = new int[m];
    int* sssRowPtrA = new int[m+1];
    int* sssColIdxA = new int[nnzA/2];
    VALUE_TYPE* sssValA = new VALUE_TYPE[nnzA/2];
    VALUE_TYPE* sssDValA = new VALUE_TYPE[m];
    
    if (csrRowPtrA[m/2] < nnzA/2)
    {
        for (int i = 0; i < m; i++)
        {
            sssDValA[i] = 0;
            sssRowPtrA[i] = num;
            for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++)
            {
                if (csrColIdxA[j] >= i)
                {
                    if (csrColIdxA[j] == i)
                        sssDValA[i] = csrValA[j];
                    break;
                }
                sssValA[num] = csrValA[j];
                sssColIdxA[num++] = csrColIdxA[j];
            }
            RowValue[i] = num - sssRowPtrA[i];
        }
    }
    else
    {
        for (int i = 0; i < m; i++)
        {
            sssDValA[i] = 0;
            sssRowPtrA[i] = num;
            for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++)
            {
                if (csrColIdxA[j] <= i)
                {
                    if (csrColIdxA[j] == i)
                        sssDValA[i] = csrValA[j];
                    continue;
                }
                sssValA[num] = csrValA[j];
                sssColIdxA[num++] = csrColIdxA[j];
            }
            RowValue[i] = num - sssRowPtrA[i];
        }        
    }
    RowValue[m-1] = num - sssRowPtrA[m-1];
    sssRowPtrA[m] = num;
    
    
    int *nRowPerRow = new int[m];
    int **Row2Row = new int*[m];

    long long conflicts = create_Conflict_Graph(sssRowPtrA, sssColIdxA, nRowPerRow, Row2Row, m);

    
    //Create DC Tree
    cout << "DCTree begin to create!\n" << flush;
    DC *DCRoot = new DC(m, partSize, numParts);
    int start_CSR = DCRoot->DC_creation(Row2Row, nRowPerRow, RowValue, m);
    cout << "DCTree created successfully!\n" << flush;
    cout << start_CSR << endl << flush;

    int *RowRev = DCRoot->DC_get_RowRev();
    int *RowPerm = DCRoot->DC_get_RowPerm();

    // Rerorder
    int* new_sssRowPtrA = new int[m+1];
    int* new_sssColIdxA = new int[nnzA/2];
    VALUE_TYPE* new_sssValA = new VALUE_TYPE[nnzA/2];
    VALUE_TYPE* new_sssDValA = new VALUE_TYPE[m];

    VALUE_TYPE* new_x = new VALUE_TYPE[m]();
    VALUE_TYPE* z = new VALUE_TYPE[m]();


    reorder_matrix(m, new_sssRowPtrA, new_sssColIdxA, new_sssValA, new_sssDValA, sssRowPtrA, sssColIdxA, sssValA, sssDValA, RowRev, RowPerm);
    reorder_vector(m, new_x, x, RowRev);

    int *new_csrRowPtrA, *new_csrColIdxA;
    VALUE_TYPE *new_csrValA;
    convert_SSS_to_CSR(m, new_sssRowPtrA, new_sssColIdxA, new_sssValA, new_sssDValA, new_csrRowPtrA, new_csrColIdxA, new_csrValA, start_CSR);

    VALUE_TYPE *y_ref = new VALUE_TYPE[m]();

    char *userArgsA[12] = {(char *)new_sssRowPtrA, (char *)new_sssColIdxA, (char *)new_sssValA, (char *)new_sssDValA, (char *)new_x, (char *)y_ref, (char *)z, (char *)new_csrRowPtrA, (char *)new_csrColIdxA, (char *)new_csrValA, (char *)&m, (char *)&start_CSR};
    char *userArgsB[12] = {(char *)new_sssRowPtrA, (char *)new_sssColIdxA, (char *)new_sssValA, (char *)new_sssDValA, (char *)y_ref, (char *)z, (char *)new_x, (char *)new_csrRowPtrA, (char *)new_csrColIdxA, (char *)new_csrValA, (char *)&m, (char *)&start_CSR};
    char *userArgsC[12] = {(char *)new_sssRowPtrA, (char *)new_sssColIdxA, (char *)new_sssValA, (char *)new_sssDValA, (char *)z, (char *)new_x, (char *)y_ref, (char *)new_csrRowPtrA, (char *)new_csrColIdxA, (char *)new_csrValA, (char *)&m, (char *)&start_CSR};
    char **userArgs[3] = { userArgsA, userArgsB, userArgsC};


    double end = omp_get_wtime();
	*pre_time = (end - start) * 1000;

    // warm up
    DCRoot->DC_traversal(SpMV_kernel_once, SpMV_kernel_twice, Final_SpMV, userArgs, K);

    *solve_time = 0;
    for(int i = 0; i < BENCH_REPEAT; i++)
    {
        // reset new_x
        reorder_vector(m, new_x, x, RowRev);
        for(int j = 0; j < m; j++) y_ref[j] = z[j] = 0;

        start = omp_get_wtime();
        DCRoot->DC_traversal(SpMV_kernel_once, SpMV_kernel_twice, Final_SpMV, userArgs, K);
        end = omp_get_wtime();
        *solve_time += (end - start) * 1000;
    }



    *solve_time /= BENCH_REPEAT;
    *gflops = 2ll * nnzA * K / (1e6 * (*solve_time));


    cout << "finishing!\n" << flush;
    
    if (K%3==0)
    {
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            y[i] = new_x[RowPerm[i]];
    }
    else if (K%3==1)
    {
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            y[i] = y_ref[RowPerm[i]];
    }
    else
    {
        #pragma omp parallel for
        for(int i = 0; i < m; i++)
            y[i] = z[RowPerm[i]];
    }

    delete[] new_x;
    delete[] y_ref;
    delete[] z;

    delete[] new_sssRowPtrA;
    delete[] new_sssColIdxA;
    delete[] new_sssValA;
    delete[] new_sssDValA;
    delete DCRoot;

    delete[] sssRowPtrA;
    delete[] sssColIdxA;
    delete[] sssDValA;
    delete[] sssValA;

    delete[] new_csrRowPtrA;
    delete[] new_csrColIdxA;
    delete[] new_csrValA;
    
    for(int i = 0; i < m; i++)
        delete[] Row2Row[i];
    delete[] Row2Row;
    delete[] nRowPerRow;
    delete[] RowValue;
}

#endif