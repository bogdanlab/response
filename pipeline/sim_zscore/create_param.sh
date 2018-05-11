grep -v start fourier_ls-chr22.bed | \
while read line
do
    region_start=$(echo $line | awk '{print $2}')
    region_stop=$(echo $line | awk '{print $3}')
    for n in 5000 50000 100000
    do
        for hsq in 0.0015 0.05 0.5
        do
            for ncau in 3
            do
                echo $n $hsq $ncau $region_start $region_stop
            done
        done
    done
done
