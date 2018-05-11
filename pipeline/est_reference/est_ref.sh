#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M shihuwen@mail
#$ -l h_data=24G,h_rt=12:00:00
#$ -j y
#$ -o ./job_out

export PYTHONPATH=/u/home/s/shihuwen/python_packages:/u/home/s/shihuwen/python_packages/lib/python2.7/site-packages
export LD_LIBRARY_PATH=/u/local/compilers/intel/current/current/lib/intel64:/u/local/compilers/intel/current/current/mkl/lib/em64t:/u/local/apps/python/2.7.3/lib:/u/local/compilers/intel/current/current/lib/intel64:/u/local/apps/hdf5/current/lib
export LD_RUN_PATH=/u/local/compilers/intel/current/current/lib/intel64
export MKLROOT=/u/local/compilers/intel/current/current/mkl
export MKL_NUM_THREADS=8
source /u/local/Modules/default/init/modules.sh
module load gcc/4.7.2
module load intel
module load python/2.7.3
mypython=/u/local/apps/python/2.7.3/bin/python

src=/u/project/pasaniuc/shihuwen/response/code/hess_test

j=$SGE_TASK_ID
line=$(head -n $j ./params.txt | tail -n 1)

n=$(echo $line | awk '{print $1}')
hsq=$(echo $line | awk '{print $2}')
ncau=$(echo $line | awk '{print $3}')
region_start=$(echo $line | awk '{print $4}')
region_stop=$(echo $line | awk '{print $5}')

for idx in $(seq 100)
do
    $mypython $src/hess.py \
        --local-hsqg /u/project/pasaniuc/shihuwen/response_result/sim_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$idx".txt --chrom 22 \
        --bfile /u/project/pasaniuc/shihuwen/posc/analysis/data/1000GP/EUR/1000GP_EUR_chr22 \
        --out /u/project/pasaniuc/shihuwen/response_result/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$idx"_step1 \
        --partition /u/project/pasaniuc/shihuwen/response/pipeline/regions/chr22_"$region_start"_"$region_stop".bed
        
    $mypython $src/hess.py \
        --prefix /u/project/pasaniuc/shihuwen/response_result/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$idx"_step1 \
        --out /u/project/pasaniuc/shihuwen/response_result/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$idx"_step2
        
    rm /u/project/pasaniuc/shihuwen/response_result/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$idx"*.gz
done
