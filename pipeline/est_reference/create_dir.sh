out_dir=/u/project/pasaniuc/shihuwen/response_result
while read line
do
    n=$(echo $line | awk '{print $1}') 
    hsq=$(echo $line | awk '{print $2}')
    ncau=$(echo $line | awk '{print $3}')
    if [ ! -d $out_dir/est_reference/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau" ]
    then
        mkdir $out_dir/est_reference/est_reference_n_"$n"_hsq_"$hsq"_ncau_"$ncau"
    fi
done < params.txt
