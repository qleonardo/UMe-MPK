#!/bin/bash
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --exclusive
#SBATCH --mem=360G
filenames=(
'audikw_1'
'bundle_adj'
'CO'
'crankseg_1'
'crankseg_2'
'dielFilterV3real'
'Ga10As10H30'
'Ga19As19H42'
'Ga3As3H12'
'Ga41As41H72'
'Ge87H76'
'Ge99H100'
'kron_g500-logn16'
'kron_g500-logn17'
'kron_g500-logn18'
'kron_g500-logn19'
'kron_g500-logn20'
'kron_g500-logn21'
'nd12k'
'nd24k'
'nlpkkt160'
'pkustk12'
'pkustk14'
'Queen_4147'
'Si41Ge41H72'
'Si87H76'
'SiO2'
)

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 5
done

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 10
done

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=28-55 --membind=1 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 15
done

