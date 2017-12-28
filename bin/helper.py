#!/usr/bin/env python
#-*- coding: utf-8 -*-

#Helper functions for TF-STAR

import sqlite3
import pandas as pd
import numpy as np
import re
import os
import sys
import random
import hashlib
import time
from collections import defaultdict
from constant import *
import subprocess
import scipy as sp
import scipy.optimize
import scipy.fftpack
import scipy.stats
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from numpy import array, empty

def check_input_id(in_file):
    '''input: gene list file
       output: gene id set and its type (ensembl, entrez, symbol)
    '''
    df = pd.read_table(in_file, sep='\n', delimiter=None, header=0, engine='python')
    ids = set(df.iloc[:,0].map(str.strip))

    ensembl_pat = re.compile(r'^ENSG\d+$')
    entrez_pat = re.compile(r'^\d+$')

    ensembl = 0
    entrez  = 0
    symbol = 0

    ids_num = len(ids) 
    for i in ids:
        if ensembl_pat.match(i):
            ensembl += 1
        elif entrez_pat.match(i):
            entrez += 1 
        else:
            symbol += 1
    max_num = max(ensembl,entrez,symbol)

    if max_num/float(ids_num) >= 0.9:
        if max_num == ensembl:
            in_type = 'ensembl'
        elif max_num == entrez:
            in_type = 'entrez'
        elif max_num == symbol:
            in_type = 'symbol'

    return ids, in_type

def background_genes_map(id_type, biotype):
    '''key: ID, value: symbol
    1   11868   14409   ENST00000456328 gene_name   +   ENSG00000223972 pseudogene
    1   30365   30503   NR_036051.1 MIR1302-2   +   100302278   noncoding
    '''
    biotype_dict = {
            '3prime_overlapping_ncrna' : 'noncoding',
            'antisense' : 'noncoding',
            'IG_C_gene' : 'protein_coding',
            'IG_C_pseudogene' : 'pseudogene',
            'IG_D_gene' : 'protein_coding',
            'IG_J_gene' : 'protein_coding',
            'IG_J_pseudogene' : 'pseudogene',
            'IG_V_gene' : 'protein_coding',
            'IG_V_pseudogene' : 'pseudogene',
            'lincRNA' : 'noncoding',
            'miRNA' : 'noncoding',
            'misc_RNA' : 'noncoding',
            'Mt_rRNA' : 'noncoding',
            'Mt_tRNA' : 'noncoding',
            'polymorphic_pseudogene' : 'protein_coding',
            'processed_transcript' : 'noncoding',
            'protein_coding' : 'protein_coding',
            'pseudogene' : 'pseudogene',
            'rRNA' : 'noncoding',
            'sense_intronic' : 'noncoding',
            'sense_overlapping' : 'noncoding',
            'snoRNA' : 'noncoding',
            'snRNA' : 'noncoding',
            'TR_C_gene' : 'protein_coding',
            'TR_D_gene' : 'protein_coding',
            'TR_J_gene' : 'protein_coding',
            'TR_J_pseudogene' : 'pseudogene',
            'TR_V_gene' : 'protein_coding',
            'TR_V_pseudogene' : 'pseudogene',
            'noncoding' : 'noncoding', # NCBI refseq gene
    }
    # biotype must in {'protein_coding', 'noncoding', 'pseudogene', 'all'}
    biotype_set = set()
    for i in biotype_dict:
        if biotype_dict[i] == biotype:
            biotype_set.add(i)
    if biotype == 'all':
        biotype_set = set(biotype_dict.keys())

    background_genes = {}

    if id_type == 'entrez':
        db_type = 'refseq'
    elif id_type == 'ensembl':
        db_type = 'ensembl'
    else:
        db_type = 'refseq'
    
    biotype_str = (', '.join('"' + item + '"' for item in biotype_set))

    conn = sqlite3.connect('./refs/references.db')
    cursor = conn.execute('SELECT geneid, symbol FROM {} WHERE biotype IN ({})'.format(db_type,biotype_str))
    for row in cursor:
        background_genes[row[0]] = row[1]

    conn.close()
    return background_genes
    

def id_converter(in_ids, id_type):
    '''return a dictionary:
    key: gene ID
    value: symbol
    '''
    out_genes = {}
    gene_to_symbol = background_genes_map(id_type, 'all')
    symbol_to_gene = {}
    for i,j in list(gene_to_symbol.items()):
        symbol_to_gene[j] = i

    for in_id in in_ids:
        if id_type == 'entrez' or id_type == 'ensembl':
            if in_id in gene_to_symbol:
                out_genes[in_id] = gene_to_symbol[in_id]
        else:
            if in_id in symbol_to_gene:
                out_genes[symbol_to_gene[in_id]] = in_id

    if len(out_genes) == 0:
        print('Please check you input, only ensembl IDs, entrez IDs and gene symbols are allowed!')
        sys.exit(1)
    return out_genes

def removefile(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise
##############################################################################
# ------------------------------------#
# P-value estimation using KDE        #
# ------------------------------------#
def _bandwidth(data, N=None, MIN=None, MAX=None):
    '''
    An implementation of the kde bandwidth selection method outlined in:

    Z. I. Botev, J. F. Grotowski, and D. P. Kroese. Kernel density
    estimation via diffusion. The Annals of Statistics, 38(5):2916-2957, 2010.

    Based on the implementation in Matlab by Zdravko Botev.

    Daniel B. Smith, PhD
    https://github.com/Daniel-B-Smith/KDE-for-SciPy/blob/master/kde.py
    Updated 1-23-2013
    '''
    # Parameters to set up the mesh on which to calculate
    N = 2**14 if N is None else int(2**sp.ceil(sp.log2(N)))
    if MIN is None or MAX is None:
        minimum = min(data)
        maximum = max(data)
        Range = maximum - minimum
        MIN = minimum - Range/10 if MIN is None else MIN
        MAX = maximum + Range/10 if MAX is None else MAX

    # Range of the data
    R = MAX-MIN

    # Histogram the data to get a crude first approximation of the density
    M = len(data)
    DataHist, bins = sp.histogram(data, bins=N, range=(MIN,MAX))
    DataHist = DataHist/M
    DCTData = scipy.fftpack.dct(DataHist, norm=None)

    I = [iN*iN for iN in range(1, N)]
    SqDCTData = (DCTData[1:]/2)**2

    # The fixed point calculation finds the bandwidth = t_star
    guess = 0.1
    try:
        t_star = scipy.optimize.brentq(__fixed_point, 0, guess, 
                                       args=(M, I, SqDCTData))
    except ValueError:
        print('Oops!')
        return None

    # Smooth the DCTransformed data using t_star
    SmDCTData = DCTData*sp.exp(-sp.arange(N)**2*sp.pi**2*t_star/2)
    # Inverse DCT to get density
    density = scipy.fftpack.idct(SmDCTData, norm=None)*N/R
    mesh = [(bins[i]+bins[i+1])/2 for i in range(N)]
    bandwidth = sp.sqrt(t_star)*R

    density = density/sp.trapz(density, mesh)

    # return bandwidth, mesh, density
    return bandwidth

def __fixed_point(t, M, I, a2):
    l=7
    I = sp.float128(I)
    M = sp.float128(M)
    a2 = sp.float128(a2)
    f = 2*sp.pi**(2*l)*sp.sum(I**l*a2*sp.exp(-I*sp.pi**2*t))
    for s in range(l, 1, -1):
        K0 = sp.prod(list(range(1, 2*s, 2)))/sp.sqrt(2*sp.pi)
        const = (1 + (1/2)**(s + 1/2))/3
        _time=(2*const*K0/M/f)**(2/(3+2*s))
        f=2*sp.pi**(2*s)*sp.sum(I**s*a2*sp.exp(-I*sp.pi**2*_time))
    return t-(2*M*sp.sqrt(sp.pi)*f)**(-2/5)

def scipy_bandwidth(data):
    kernel = scipy.stats.gaussian_kde(data)
    f = kernel.covariance_factor()
    bw = f * data.std()
    return bw


def kde(data, bw):
    #instantiate and fit the KDE model
    kde = KernelDensity(bandwidth=bw,kernel='gaussian')
    kde.fit(data[:,np.newaxis])
    
    # obtain the range of data 
    # and define the range of distribution
    minimum = min(data)
    maximum = max(data)
    Range = maximum - minimum
    MIN = round(minimum - Range/2) 
    MAX = round(maximum + Range/2)

    data_points = np.around(np.arange(MIN, MAX, 0.01, dtype=np.float64),3)

    # the log of the probability density
    logprob = kde.score_samples(data_points[:,None])
    prob = np.exp(logprob)

    '''
    plt.fill_between(data_points, np.exp(logprob), alpha=0.5)
    plt.plot(data, np.full_like(data, -0.01), '|k', markeredgewidth=1)
    plt.show()
    '''
    #print(sp.trapz(np.exp(logprob), data_plot))
    return prob, data_points # y, x

def pvalue_calc(input_value, array_x, array_y):
    '''right-side one-tailed pvalue
    '''
    minimum = min(array_x)
    maximum = max(array_x)
    if input_value <= minimum:
        return 1.0
    if input_value > maximum:
        return 0.0
    _x = np.array(array_x)
    #_x_str = np.array(map(str,array_x))
    _y = np.array(array_y)
    density = dict(list(zip(_x, _y)))
    x_values = _x[_x >= input_value]
    y_values = [density[i] for i in x_values]
    return sp.trapz(y_values, x_values)

def best_bandwidth_by_cv(data):
    min_bound = np.mean(data)/10.0 + 1
    max_bound = np.mean(data)/10.0 - 1
    bandwidths = np.arange(min_bound,max_bound,0.1)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                       {'bandwidth':bandwidths},
                       cv=10,
                       n_jobs=10)
    grid.fit(data[:,None])
    return grid.best_params_

def correct_pvalues_for_multiple_testing(pvalues, correction_type="Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    pvalues = array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


def hypergeometric_pvalue(k, M, n, N):
    ''' M: total number of genes (population size) 
        N: the number of input genes (sample size)
        n: total number of genes binding by the test TF (sucesses size in the population) 
        k: the number of input genes binding by the test TF (sucesses size in the sample)  
    # this is equivalent to "1 - phyper(k-1, n, M-n, N)"
    '''
    p = scipy.stats.hypergeom.sf(k-1, M, n, N)
    return '{0:.5g}'.format(p)



if __name__ == '__main__':
    #print(id_converter({'IKZF1','STAT4'}, 'symbol'))
    #a,b = refseq_map('entrez', 'protein_coding', '2000:200', '200:2000')
    #print(overlap_peaks(a, 'GM12878.encode.bed'))
    data = np.arange(1, 100, 0.1)
    bw = _bandwidth(data)
    print(bw)
    bw = scipy_bandwidth(data)
    print(bw)
    #bw = best_bandwidth_by_cv(data)
    #print(bw)
    #print(kde(data, bw))
