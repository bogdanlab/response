export PYTHONPATH=/u/home/s/shihuwen/.local/lib/python3.6/site-packages:/u/project/pasaniuc/shihuwen/software/ldetect/lib/python3.6/site-packages:/u/project/pasaniuc/shihuwen/software/commanderline-0.2/lib/python3.6/site-packages
export LD_LIBRARY_PATH=/u/local/apps/python/3.1.2/lib:/u/local/apps/python/3.6.1/lib

python3 plot_sim_result.py \
    --result ../sim_results/insample_all_eig.txt --hsq 0.0015 \
    --out insample_0.0015

python3 plot_sim_result.py \
    --result ../sim_results/insample_all_eig.txt --hsq 0.05 \
    --out insample_0.05

python3 plot_sim_result.py \
    --result ../sim_results/insample_all_eig.txt --hsq 0.5 \
    --out insample_0.5
