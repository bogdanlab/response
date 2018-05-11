import numpy as np
import pandas as pd
import argparse, math, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec


matplotlib.rc('font', family='sans-serif') 
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['ytick.labelsize'] = 18
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# main function
def main():
    
    args = get_command_line()

    result = pd.read_table(args.result, delim_whitespace=True)

    region_stat = result[['P', 'START', 'STOP']].drop_duplicates()
    avg_nsnp = int(np.mean(region_stat['P']))
    avg_width = int(np.mean(region_stat['STOP']-region_stat['START']))

    fig, ax = plt.subplots(1,1,figsize=(7.5, 7.5)) 
 
    ax.boxplot(result['EST'][(result['N']==50000)&(result['HSQ']==args.hsq)].tolist(),
               positions=[1], widths=0.5, showmeans=True)

    ax.boxplot(result['EST'][(result['N']==75000)&(result['HSQ']==args.hsq)].tolist(),
               positions=[2], widths=0.5, showmeans=True)

    ax.boxplot(result['EST'][(result['N']==100000)&(result['HSQ']==args.hsq)].tolist(),
               positions=[3], widths=0.5, showmeans=True)

    ax.set_xlim([0.5, 3.5])
    ax.plot([0, 4], [args.hsq, args.hsq], 'k--', 
        label=r'simulated $\boldsymbol{\beta}^T\boldsymbol{V}\boldsymbol{\beta}$ ' + '({})'.format(args.hsq))
    ax.set_xticks([1,2,3])
    ax.set_xticklabels([50000, 75000, 100000])
    ax.set_xlabel('sample size', fontsize=20)
    ax.set_ylabel(r'estimated $\boldsymbol{\beta}^T\boldsymbol{V}\boldsymbol{\beta}$', fontsize=20)
    ax.set_title('using in-sample LD \n average number of SNPs: {} \n average region size: {}'.format(avg_nsnp, avg_width), fontsize=18)

    plt.legend(loc=1, fontsize=18, frameon=False)

    plt.savefig('{}.pdf'.format(args.out_file), bbox_inches='tight')

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='plot convergence')

    parser.add_argument('--result', dest='result', type=str,
        help='result file', required=True)

    parser.add_argument('--hsq', dest='hsq', type=float,
        help='heritability', required=True)

    parser.add_argument('--out', dest='out_file', type=str,
        help='output file name', required=True)
   
    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
