#!/bin/bash
cd ../../../matrix/

for matrix in *
do
    cd ../MPK/DC_MPK_1.12/MPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 5
    cd ../../../matrix/
done

for matrix in *
do 
    cd ../MPK/DC_MPK_1.12/MPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 10
    cd ../../../matrix/
done

for matrix in *
do 
    cd ../MPK/DC_MPK_1.12/MPK/
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 15
    cd ../../../matrix/
done
