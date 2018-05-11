njobs=$(cat params.txt | wc  -l)
qsub -t 1-$njobs:1 est_1kg.sh
