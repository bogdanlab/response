njobs=$(cat params.txt | wc  -l)
qsub -t 1-$njobs:3 sim_zsc.sh
