#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M shihuwen@mail
#$ -l h_data=12G,h_rt=6:00:00
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

src=/u/project/pasaniuc/shihuwen/response/code/sim
params=/u/project/pasaniuc/shihuwen/response/pipeline/sim_zscore/params.txt

i=$SGE_TASK_ID
for j in $(seq $i $(($i+2)))
do
    line=$(head -n $j $params | tail -n 1)

    n=$(echo $line | awk '{print $1}')
    hsq=$(echo $line | awk '{print $2}')
    ncau=$(echo $line | awk '{print $3}')
    region_start=$(echo $line | awk '{print $4}')
    region_stop=$(echo $line | awk '{print $5}')

    echo $n $hsq $ncau $region_start $region_stop

    $mypython $src/sim_gwas.py \
        --n $n --hsq $hsq --num_sim 50 --ncau $ncau \
        --region $region_start $region_stop \
        --legend /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb_chr22_maf/22.bim \
        --bfile /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb_chr22_maf/22 \
        --out /u/project/pasaniuc/shihuwen/response_result/sim_zsc/sim_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_
done 
