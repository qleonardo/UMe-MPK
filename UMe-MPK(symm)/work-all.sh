#!/bin/bash
cd ../../../matrix/

for matrix in *
do
    cd ../MPK/SymDC_MPK_1.12/SymMPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 5
    cd ../../../matrix/
done

for matrix in *
do 
    cd ../MPK/SymDC_MPK_1.12/SymMPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 10
    cd ../../../matrix/
done

for matrix in *
do 
    cd ../MPK/SymDC_MPK_1.12/SymMPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 15
    cd ../../../matrix/
done
