njobs=$(cat params.txt | wc  -l)
qsub -t 1-$njobs:1 est_uk10k.sh
