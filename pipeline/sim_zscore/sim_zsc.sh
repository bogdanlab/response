src=/u/project/pasaniuc/shihuwen/response/code/sim
params=/u/project/pasaniuc/shihuwen/response/pipeline/sim_zscore/params.txt

while read line
do
    n=$(echo $line | awk '{print $1}')
    hsq=$(echo $line | awk '{print $2}')
    ncau=$(echo $line | awk '{print $3}')
    region_start=$(echo $line | awk '{print $4}')
    region_stop=$(echo $line | awk '{print $5}')

    echo $n $hsq $ncau $region_start $region_stop

    python $src/sim_gwas_pop_chrom_ukb_cvg.py \
        --n $n --hsq $hsq --num_sim 50 --ncau $ncau \
        --region $region_start $region_stop \
        --legend /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb_chr22_maf/22.bim \
        --bfile /u/project/pasaniuc/shihuwen/posc/analysis/data/ukb_chr22_maf/22 \
        --out /u/project/pasaniuc/shihuwen/response_result/sim_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_
done < $params
