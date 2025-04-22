#!/bin/bash
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --exclusive
#SBATCH --mem=360G
cd ../../../matrix/
filenames=(
'af_shell10'
# 'audikw_1'
# 'cant'
# 'G3_circuit'
# 'Hook_1498'
# 'inline_1'
# 'ldoor'
# 'nlpkkt120'
# 'pwtk'
# 'Serena'
# 'shipsec1'

# 'cfd2'
# 'parabolic_fem'
# 'offshore'
# 'bmw7st_1'
# 'thermal2'
# 'gearbox'
# 'crankseg_1'
# 'F1'
# 'Fault_639'
# 'Geo_1438'
# 'bone010'
# 'dielFilterV3real'
)

# for matrix in ${filenames[@]}
# do 
#     cd ../qhz/MPK/MPK/
#     OMP_NUM_THREADS=24 numactl --physcpubind=0-23 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 1
#     cd ../../../matrix/
# done

for matrix in ${filenames[@]}
do 
    cd ../qhz/MPK/MPK/
    OMP_NUM_THREADS=24 numactl --physcpubind=0-23 --membind=0 ./MAIN ../../../matrix/$matrix/$matrix.mtx $matrix 5
    cd ../../../matrix/
done

# for matrix in ${filenames[@]}
# do 
#     cd ../qhz/MPK/MPK/
#     OMP_NUM_THREADS=24 numactl --physcpubind=0-23 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 10
#     cd ../../../matrix/
# done