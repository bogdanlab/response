#!/usr/bin/python

from optparse import OptionParser

import pandas as pd
import numpy as np
import numpy.linalg as nplinalg
import random, math, sys, argparse
from pysnptools.snpreader import Bed
import os.path

# main function
def main():
    
    # get command line
    args = get_command_line()

    # load legend and genotype data
    legend = pd.read_table(args.legend, sep='\t', header=None)
    legend.columns = ['CHR', 'SNP', 'CM', 'BP', 'A2', 'A1']

    # find the index of the genotypes
    incl_snps_idx = legend[(legend['BP']>=args.region[0]) & \
                           (legend['BP']<args.region[1])]

    print(incl_snps_idx)

    sys.exit(0)

    # load genotype data
    genotype = Bed(args.bfile, count_A1=False)
    genotype = genotype.read().val

    # impute missing data using sample average
    nanidx = np.where(np.isnan(genotype))
    mean_geno = np.nanmean(genotype, axis=0)
    genotype[nanidx] = mean_geno[nanidx[1]] 
    genotype = genotype1[0:args.n, ]

    # standardize genotypes
    geno_std = std_geno(genotype)

    # get the ld matrix
    num_snps = legend.shape[0]

    # iterate through total number of simulations
    for i in range(args.num_sim):

        # draw causal snps
        cau = draw_causals(args.n, num_snps)

        # simulate phenotype
        phe,beta,phe_g,phe_e = sim_pheno(geno_std, cau, args.hsq, num_snps)

        # simulate z-scores
        zsc = sim_zsc(geno_std, phe)

        # write out result
        write_zsc(legend, zsc, args.n1, '%s%d.pop1'.format(args.out, i+1))

# standardize genotypes
def std_geno(geno):
    geno -= geno.mean(axis=0)
    geno /= geno.std(axis=0)
    return np.nan_to_num(geno)

# draw causal snps
def draw_causals(ncau, nsnps):
    
    all_index = np.array(range(nsnps))
    cau_idx = np.random.choice(range(nsnps), size=ncau, replace=False)
    cau = all_index[cau_idx]

    return cau

# simulate phenotype
def sim_pheno(geno, cau, hsq, nsnp):
    
    # draw effect sizes for causal snps
    geno_cau = geno[:,cau]
    num_cau = geno_cau.shape[1]
    beta = np.random.normal(loc=0.0, scale=np.sqrt(hsq/num_cau), size=num_cau)
    beta_full = np.zeros(nsnp) 
    beta_full[cau] = beta

    # simulate the genetic component
    phe_g = np.dot(geno_cau, beta)
    var_phe_g = np.var(phe_g)
    
    ratio = hsq / var_phe_g
    phe_g = phe_g * np.sqrt(ratio)

    # add the environment component
    num_indv = geno_cau.shape[0]
    phe_e = np.random.normal(loc=0.0, scale=np.sqrt(1.0-hsq), size=num_indv)
    var_phe_e = np.var(phe_e)
    ratio = (1-hsq) / var_phe_e
    phe_e = phe_e * np.sqrt(ratio)

    # simulate the trait
    phe = phe_g + phe_e

    return phe,beta_full,phe_g,phe_e

# simulate z scores
def sim_zsc(geno, phe):
    num_indv = geno.shape[0]
    zsc = np.dot(geno.T, phe)/np.sqrt(num_indv)
    return zsc

# write out result
def write_zsc(legend, zsc, beta, cau, n, file_nm):
    cau_set = set(np.ndarray.tolist(cau))
    file_nm = open(file_nm, 'w')
    file_nm.write('SNP\tCHR\tBP\tA2\tA1\tZ\tN\n')
    for i in range(zsc.shape[0]):
        status = 0
        if i in cau_set:
            status = 1
        file_nm.write('%s\t%d\t%d\t%s\t%s\t%f\t%d\n' % (
            legend['SNP'][i], legend['CHR'][i], legend['BP'][i],
            legend['A2'][i], legend['A1'][i], zsc[i], n))
    file_nm.close()

# get command line
def get_command_line():
    parser = argparse.ArgumentParser(description='simulate summary stats')

    parser.add_argument('--legend', dest='legend', type=str,
        help='legend file', required=True)
    parser.add_argument('--bfile', dest='bfile', type=str,
        help='genotype file', required=True)

    parser.add_argument('--ncau', dest='ncau', type=float,
        help='number of causal variants', required=True)

    parser.add_argument('--hsq', dest='hsq', type=float,
        help='heritability of the trait', required=True)

    parser.add_argument('--n', dest='n', type=float,
        help='sample size of the gwas', required=True)

    parser.add_argument('--region', dest='region', type=int, nargs=2,
        help='start and stop position of the region', required=True)

    parser.add_argument('--snps', dest='snps', type=str,
        help='a list of snps to include', required=True)

    parser.add_argument('--num_sim', dest='num_sim', type=int,
        help='number of replications to generate', required=True)

    parser.add_argument('--out', dest='out', type=str,
        help='output', required=True)
    
    args = parser.parse_args()
    
    return args

# run main
if __name__ == "__main__":
    main()
