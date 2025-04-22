#!/bin/bash
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --exclusive
#SBATCH --mem=360G
filenames=(
#'ASIC_320k'
#'as-Skitter'
#'audikw_1'
#'bundle_adj'
#'com-LiveJournal'
#'com-Orkut'
#'com-Youtube'
#'coPapersCiteseer'
#'crankseg_2'
#'dielFilterV3real'
#'FullChip'
#'Ga10As10H30'
#'Ga19As19H42'
#'Ga3As3H12'
#'Ga41As41H72'
#'GaAsH6'
#'Ge87H76'
#'Ge99H100'
#'gupta2'
#'gupta3'
#'hollywood-2009'
#'ins2'
#'kron_g500-logn16'
#'kron_g500-logn17'
#'kron_g500-logn18'
#'kron_g500-logn19'
#'kron_g500-logn20'
#'kron_g500-logn21'
#'lp1'
#'nd12k'
#'nd24k'
#'nlpkkt120'
#'nlpkkt160'
#'Queen_4147'
#'rajat30'
#'Si41Ge41H72'
#'Si87H76'
#'SiO2'
'cfd2'
)

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 5
done

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 10
done

for matrix in ${filenames[@]}
do 
    OMP_NUM_THREADS=28 numactl --physcpubind=0-27 --membind=0 ./main ../../../matrix/$matrix/$matrix.mtx $matrix 15
done
