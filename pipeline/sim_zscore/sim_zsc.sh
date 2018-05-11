export PYTHONPATH=/u/home/s/shihuwen/.local/lib/python3.6/site-packages:/u/project/pasaniuc/shihuwen/software/ldetect/lib/python3.6/site-packages:/u/project/pasaniuc/shihuwen/software/commanderline-0.2/lib/python3.6/site-packages
export LD_LIBRARY_PATH=/u/local/apps/python/3.1.2/lib:/u/local/apps/python/3.6.1/lib

src=/u/project/pasaniuc/shihuwen/response/code/sim

n=5000
hsq=0.05
ncau=3

python3 $src/sim_gwas_pop_chrom_ukb_cvg.py \
    --n $n --hsq $hsq --num_sim 100 \
    --legend /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb/22.bim \
    --bfile /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb/22 \
    --out $out_dir/sim_gwas

