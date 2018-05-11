result_dir=/u/project/pasaniuc/shihuwen/response_result/est_insample
echo N HSQ P NCAU START STOP EST

while read line
do
    n=$(echo $line | awk '{print $1}')
    hsq=$(echo $line | awk '{print $2}')
    ncau=$(echo $line | awk '{print $3}')
    region_start=$(echo $line | awk '{print $4}')
    region_stop=$(echo $line | awk '{print $5}')

    for i in $(seq 50)
    do
        out=$result_dir/est_insample_n_"$n"_hsq_"$hsq"_ncau_"$ncau"/sim_gwas_"$region_start"_"$region_stop"_"$i"_step2.log
        nsnp=$(grep Using $out | awk '{print $3}')
        est=$(grep Total $out | awk '{print $5}')
        echo $n $hsq $nsnp $ncau $region_start $region_stop $est
    done
    
done < params.txt
