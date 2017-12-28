#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# __author__ = 'T-Y Wang'
# ------------------------------------
# Python Module
# ------------------------------------

import os
import sys
import argparse
import subprocess
import multiprocessing
import random
import tqdm
import datetime
import logging
import errno
import re
import json
import io
import openpyxl
import numpy as np
from time import sleep,time
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import hashlib
import glob
import itertools
import scipy as sp
import scipy.optimize
import scipy.fftpack
import scipy.stats
import time
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import sqlite3
#from helper import check_input_id, background_genes_map, id_converter, removefile   
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)
warnings.simplefilter(action = "ignore", category = UserWarning)

# --------------------------------------------- #
# define constants                              #
# --------------------------------------------- #
MAX_WORKERS = multiprocessing.cpu_count()

CHRS =  {'1'  : 249250621,
         '10' : 135534747,
         '11' : 135006516,
         '12' : 133851895,
         '13' : 115169878,
         '14' : 107349540,
         '15' : 102531392,
         '16' : 90354753,
         '17' : 81195210,
         '18' : 78077248,
         '19' : 59128983,
         '2'  : 243199373,
         '20' : 63025520,
         '21' : 48129895,
         '22' : 51304566,
         '3'  : 198022430,
         '4'  : 191154276,
         '5'  : 180915260,
         '6'  : 171115067,
         '7'  : 159138663,
         '8'  : 146364022,
         '9'  : 141213431,
         'MT' : 16569,
         'X'  : 155270560,
         'Y'  : 59373566,
    }

io.DEFAULT_BUFFER_SIZE = 0 

# ------------------------------------ #
# helper functions                       #
# ------------------------------------ #

from helper import check_input_id, background_genes_map, id_converter, removefile
from helper import _bandwidth, __fixed_point, kde, pvalue_calc, best_bandwidth_by_cv, scipy_bandwidth  
from helper import hypergeometric_pvalue, correct_pvalues_for_multiple_testing

###################################################################################
# -------------------------------------------------------#
# Main functions
# -------------------------------------------------------# 
def reference_map(id_type, biotype, upstream, downstream):
    '''input: id_type of input genes and selected biotype  
       output: reference gene file for TF peaks overlapping 
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

    if id_type == 'entrez':
        db_type = 'refseq'
    elif id_type == 'ensembl':
        db_type = 'ensembl'
    else:
        db_type = 'refseq'
    
    # generate unique ID for upstream_downstream_tss.bed and overlap.bed files
    # the two bed files use the identical ID 
    random.seed(time.time())
    iid = hashlib.sha224(str(random.random()).encode()).hexdigest()
    """add complicated upstream and downstream definition
       upstream: 2000:500; downstream : 100:2000
       upstream: 2000; downstream: 2000 == upstream: 2000:0; downstream: 2000:0
    """
    # checking input parameters: upstream and downstream
    upstream = str(upstream)
    downstream = str(downstream)
    pat = re.compile(r'^0:|:0$')
    if upstream.count(':') == 1:
        mat = pat.search(upstream) 
        if mat:
            sys.exit(1)
    if downstream.count(':') == 1:
        mat = pat.search(downstream) 
        if mat:
            sys.exit(1)

    tmp_folder = './tmp/'
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
    tss_bed_filename = tmp_folder + 'tss_upstream_' + str(upstream.replace(':','_')) + '_downstream_' + str(downstream.replace(':','_')) + '.' +iid + '.bed'

    tss_bed = io.open(tss_bed_filename, 'w', encoding='utf-8')
    stored_region_dict = {}
    gene_tss_number = defaultdict(int)
    
    # connect sqlite3 reference database
    biotype_str = (', '.join('"' + item + '"' for item in biotype_set))

    conn = sqlite3.connect('./refs/references.db')
    cursor = conn.execute('SELECT chr,start,end,transcript,strand,geneid FROM {} WHERE biotype IN ({})'.format(db_type,biotype_str))
    for row in cursor:
        chrm  = row[0]
        start = row[1]
        end   = row[2]
        transcript = row[3]
        strand = row[4]
        gene = row[5]

        #  complicated upstream and downstream definition
        if ':' in upstream and ':' in downstream:
            ups = list(map(int, upstream.split(':')))
            downs = list(map(int, downstream.split(':')))
            ups.sort()
            downs.sort()
            up_s, up_b = ups
            down_s, down_b = downs
                
            gene_tss_number[gene] = gene_tss_number[gene] + 1
            if strand == '+':
                # It's certain that j < k and l < m
                # ---j----k---TSS---l---m---
                j = int(start) - int(up_b)
                k = int(start) - int(up_s)
                if not j == k:
                    if j <= 0 and k <= 0:
                        pass
                    elif j <=0 and k > 0:
                        new_start = 0
                        new_end = k
                    elif j > 0 and k > 0:
                        new_start = j
                        new_end = k

                m = int(start) + int(down_b)
                l = int(start) + int(down_s)
                if not m == l:
                    if m >= CHRS[chrm] and l >= CHRS[chrm]:
                        pass 
                    elif m >= CHRS[chrm] and l < CHRS[chrm]:
                        new_start = l
                        new_end = CHRS[chrm]
                    elif m < CHRS[chrm] and l < CHRS[chrm]:
                        new_start = l
                        new_end = m

            elif strand == '-':
                # It's certain that j < k and l < m
                # ---l---m---TSS---j---k--- 
                k = int(end) + int(up_b)
                j = int(end) + int(up_s)
                if not j == k:
                    if j >= CHRS[chrm] and k >= CHRS[chrm]:
                        pass 
                    elif k >= CHRS[chrm] and j < CHRS[chrm]:
                        new_start = j
                        new_end = CHRS[chrm]
                    elif j < CHRS[chrm] and k < CHRS[chrm]:
                        new_start = j
                        new_end = k


                l = int(end) - int(down_b)
                m = int(end) - int(down_s)
                if not l == m:
                    if l <= 0 and m <= 0:
                        pass
                    elif l <=0 and m > 0:
                        new_start = 0
                        new_end = m
                    elif l > 0 and m > 0:
                        new_start = l
                        new_end = m

            stored_region = chrm + ':' + strand + ':' + str(new_start) + ':' + str(new_end)
            if not stored_region in stored_region_dict:
                stored_region_dict[stored_region] = True
                tss_bed.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrm,new_start,new_end,transcript,'.',strand,gene))
        # simple upstream and downstream definition       
        else:
            gene_tss_number[gene] = gene_tss_number[gene] + 1
            if strand == '+':
                new_start = 0 if int(start) - int(upstream) <=0 else int(start) - int(upstream)
                new_end   = CHRS[chrm] if int(start) + int(downstream) >= CHRS[chrm] else int(start) + int(downstream)
            elif strand == '-':
                new_start = 0 if int(end) - int(downstream) <= 0 else int(end) - int(downstream)
                new_end   = CHRS[chrm] if int(end) + int(upstream) >= CHRS[chrm] else int(end) + int(upstream)
                
            stored_region = chrm + ':' + strand + ':' + str(new_start) + ':' + str(new_end)
            if not stored_region in stored_region_dict:
                stored_region_dict[stored_region] = True
                tss_bed.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrm,new_start,new_end,transcript,'.',strand,gene))

    tss_bed.close()        
    stored_region_dict = {} 
    conn.close()

    return tss_bed_filename, gene_tss_number

def overlap_peaks(tss_ref, tf_peak):
    # intersect TSS bed to TF ChIP-seq bed
    # dict of geneid : {TFs} 
    # TF records: TF1,TF2,TF3 ; TF1 ; TF2 ; TF4
    # gene: TF1,TF2,TF3:TF1:TF2:TF4

    tf_id = os.path.basename(tf_peak).split('.')[0]

    tmp_folder = './tmp/'
    overlap_bed = '{}.{}.overlap.bed'.format(os.path.splitext(tss_ref)[0], tf_id)
    cmd = 'bedtools intersect -wo -a {} -b {} > {}'.format(tss_ref, tf_peak, overlap_bed)
    #print(cmd)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    if error_msg:
        print(('*** Error for generating overlapping bed file: \n{}\n{}'.format(err, error_msg)))

    # overlapping bed filtering 
    #1   860192  861992  XM_005244723.1  .   +   148398  1   860180  860330  KDM5B   0   KDM5B(K562) 138
    #1   860192  861992  XM_005244723.1  .   +   148398  1   860605  860755  KDM5B   0   KDM5B(K562) 150
    gene_tf_map = {}
    with io.open(overlap_bed, 'r', encoding='utf-8') as f:
        for line in f:
            l = line.rstrip().split('\t')
            geneID     = l[6]
            tfs        = l[10]
            tss_start  = int(l[1])
            tss_end    = int(l[2])
            peak_start = int(l[8])
            peak_end   = int(l[9])
            overlap    = int(l[-1])
            flag = False
            # guarantee always including TF peaks in overlapping regions
            if tss_end - tss_start >= peak_end - peak_start:
                if overlap == peak_end - peak_start:
                    flag = True
            else:
                if overlap == tss_end - tss_start:
                    flag = True
            if flag:
                if not geneID in gene_tf_map:
                    gene_tf_map[geneID] = tfs
                else:
                    gene_tf_map[geneID] = gene_tf_map[geneID] + ':' + tfs

    #removefile(tss_ref)
    #removefile(overlap_bed) 

    print("gene TF map loaded!")
    # modified Nov. 27, 2017
    return gene_tf_map, overlap_bed

def tf_calc(tgt_genes, gene_tf_map, gene_tss_number):
    '''
    Input ----------------------------------
    tgt_genes: input gene ids
    gene_tss_num: gene to TSS number dict
    gene_tf_map: gene to tfs dict
    Output ---------------------------------

    '''
    # extract TF bindings peaks for input gene ids
    # gene => TF1:TF2:TF1,TF2:TF3
    tfs = {}
    for gene in tgt_genes:
        try:
            tfs[gene] = gene_tf_map[gene]
        except KeyError:
            pass
    # no genes found have TF peaks
    if len(tfs) == 0:
        print('No TFs found for input genes!')
        #sys.exit(1)
        #_tfs = geneID_tf_map.values()
        #if len(_tfs) == 1:
        #    return {_tfs[0]:0}
        #else:
        return None
    # input genes have TF binding peaks
    else:
        '''
        gene => TF1:TF2:TF1,TF2:TF3
        [(TF1,),(TF2,),(TF1,TF2,),(TF3,)]
        '''
        count_total = Counter()
        for gene in tfs:
            # transcript number of the gene
            num = gene_tss_number[gene]
            tf_list = []
            peaks = tfs[gene].split(':')
            for peak in peaks:
                l = peak.split(',')
                l.sort()
                tf_list.append(tuple(l))

            count = Counter()
            # tf_list: [(TF1,),(TF2,),(TF1,TF2,),(TF3,)]
            for tf in tf_list:
                count.update(tf)
            for i in count:
                count[i] = count[i]/float(num)
            # update TFs for the input gene list
            count_total.update(count)
        # Counter({'TF1': 1.25, 'TF2': 0.416, 'TF3': 0.416})
        return count_total

def aggregate_tf_count(observed, distribution, method='empirical'):
    '''P-values calculation based on observation and permutated distribution
       available methods: "empirical", "kde"
    '''
    tf_names = list(observed.keys())
    permutation_count = dict((k, 0) for k in tf_names)
    permutation_p = dict((k, 0) for k in tf_names)
    # TF counts in distribtion
    tf_numbers = dict((k, []) for k in tf_names)
    
    tf_numbers_mean = dict((k, 0) for k in tf_names)
    tf_numbers_std = dict((k, 0) for k in tf_names)
    fold_change = dict((k, 0) for k in tf_names)

    for permutation in distribution:
        # for every permutation, store TFs and their counts in a dict
        tf_permutation = {}
        for i in tf_names:
            if i in permutation:
                tf_permutation.update({i: permutation[i]})
            else:
                tf_permutation.update({i: 0.0})

        for tf in tf_permutation:
            tf_numbers[tf].append(tf_permutation[tf])
            if observed[tf] <= tf_permutation[tf]:
                permutation_count[tf] = permutation_count[tf] + 1 

    permutation_result = {}

    for i in tf_names:
        if method == 'empirical': 
            permutation_p[i] = float(permutation_count[i])/len(distribution)
        elif method == 'kde':
            array_of_tf_numbers = np.asarray(tf_numbers[i])
            bw = _bandwidth(array_of_tf_numbers)
            #bw = best_bandwidth_by_cv(array_of_tf_numbers) 
            zeros = len(array_of_tf_numbers[np.where(array_of_tf_numbers==0.0)])
            y, x = kde(data=array_of_tf_numbers, bw=bw)
            overall = sp.trapz(y,x)
            __p = pvalue_calc(observed[i], x, y)/overall

            # print zeros, i
            # detect sample distribution or just remove the zeros?
            if zeros >= 20 or __p > 1.0:
                permutation_p[i] = float(permutation_count[i])/len(distribution)
            else:
                permutation_p[i] = __p
        else:
            pass
        tf_numbers_mean[i] = round(np.mean(tf_numbers[i]),4)
        if np.sum(tf_numbers[i]) != 0.0:
            tf_numbers_std[i] = round(np.std(tf_numbers[i]), 4)
            fold_change[i] = round(observed[i]/float(np.mean(tf_numbers[i])), 4)
        else:
            tf_numbers_std[i] = 'NaN'
            fold_change[i] = 'NaN'
        # observed_num, distribution_mean, SD, fold_change, pvalue 
        # Nov. 28, 2017 modified the order of output
        permutation_result[i] = '{:.5g}\t{}\t{}\t{}\t{}'.format(observed[i], tf_numbers_mean[i], tf_numbers_std[i], fold_change[i], permutation_p[i])

    return permutation_result


def permutation(background_genes, gene_tf_map, gene_tss_number, tgt_genes, n, threads, log_file, method):
    '''P-value calculation via permutation procedure
    '''
    print("Permutation begins!")
    print((datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))

    threads = 1

    #workers = min(MAX_WORKERS, threads)
    workers = 1
    #background_genes = gene_map.keys()
    observed_tf_count = tf_calc(tgt_genes, gene_tf_map, gene_tss_number)
    if observed_tf_count == None:
        #'observed_count\tpermutation_avg_count\tpermutation_count_std\tfold_change\tpermutation_pvalue'
        return {'placeholder':'0\tNA\tNA\tNA\t1.0'}

    tf_count_distribution = []

    count = 0
    if log_file:
        log = io.open(log_file, 'w', encoding='utf-8')

    with ThreadPoolExecutor(max_workers=workers) as executor:
        tf_stats = {}
        for i in range(n):
            permut_genes = random.sample(background_genes, len(tgt_genes))
            future = executor.submit(tf_calc, permut_genes, gene_tf_map, gene_tss_number)
            tf_stats[future] = i

        # aggregate the permutation results
        for future in tqdm.tqdm(as_completed(tf_stats),total=n):
            count = count + 1

            if log_file:
                percent = '{}'.format(int(round(100*count/float(n))))
                message = '{} permutation(s) processed.'.format(count)
                data = {"percent": int(percent), "message": message}
                #log.flush()
                log.seek(0)
                # output to JSON
                log.write(str(json.dumps(data, ensure_ascii=False)))
            try:
                res = future.result()
                if res:
                    tf_count_distribution.append(res)
            except Exception as exc:
                error_msg = 'error: ' + str(exc)
            else:
                error_msg = ''
            if error_msg:
                idx = tf_stats[future]
                print(('*** Error for {}: {}'.format(idx, error_msg)))
    if log_file:
        log.close()

    permutation_result = aggregate_tf_count(observed_tf_count, tf_count_distribution, method)

    print("Permutation accomplished!")
    print((datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))

    return permutation_result


def output_result(results_list,output_prefix,output_suffix,pairwise=True,biosample='',sampling=True): 
    '''
    Non-pairwise permutation:
        pairwise=False, biosample='GM12878'
    Pairwise permutation:
        pairwise=True, biosample=''
    '''
    if sampling:
        header = 'TF\tBiosample\tobserved_count\tpermutation_avg_count\tpermutation_count_std\tfold_change\tpermutation_pvalue\tFDR\n'
    else:
        header = 'TF\tBiosample\tbinding_NO_input_genes\tbinding_NO_total_genes\thypergeometric_pvalue\tFDR\n'
    if output_suffix == 'csv':
        header = header.replace('\t',',')
    output_file = '{}.{}'.format(output_prefix, output_suffix)

    tmp_list = [] 
    if not pairwise:
        for tf in results_list:
            tmp_list.append('{}\t{}\t{}'.format(tf,biosample,results_list[tf]))
    else:
        tmp_list = results_list

    pvalues = [ float(it.split('\t')[-1]) for it in tmp_list ] 
    fdrs = correct_pvalues_for_multiple_testing(pvalues)
    
    output_list = [] 
    for i,j in zip(tmp_list, fdrs):
        output_list.append('{0}\t{1:.5g}'.format(i,j))

    with io.open(output_file, 'w', encoding="utf-8") as f:
        f.write(header)    
        for item in output_list:
            if output_suffix == 'csv':
                item = item.replace('\t', ',')
            f.write('{}\n'.format(item))

    return output_file

def pairwise_permutation(tss_ref, tf_peak, background_genes, gene_tss_number, tgt_genes, n, method):
    ''' A simple combination function of overlap_peaks and permutation
        for pairwise TFs permutation
    '''
    gene_tf_map, overlap_bed = overlap_peaks(tss_ref, tf_peak)
    threads = 1
    log_file = False
    res = permutation(background_genes, gene_tf_map, gene_tss_number, tgt_genes, n, threads, log_file, method)
    #removefile(tss_ref)
    #removefile(overlap_bed)
    return res

def hypergeometric_calculation(tss_ref, tf_peak, background_genes, tgt_genes, tf):
    '''
    calculate hypergeometric pvalue for individual TF peak
    obtain empty string if no binding TFs for input genes, 
    binding_NO_input_genes, binding_NO_total_genes, hypergeometric_pvalue
    '''
    dict_k = defaultdict(int)
    dict_n = defaultdict(int)
    gene_tf_map, overlap_bed = overlap_peaks(tss_ref, tf_peak)
    dict_n.update(gene_binding_count4hypergeometric(in_genes=background_genes, gene_tf_map=gene_tf_map))
    dict_k.update(gene_binding_count4hypergeometric(in_genes=tgt_genes, gene_tf_map=gene_tf_map))

    M = len(background_genes)
    N = len(tgt_genes)
    
    # no hit
    if len(dict_k.keys()) == 0:
        pvalue = 1.0
        return ''
    pvalue = hypergeometric_pvalue(dict_k[tf], M, dict_n[tf], N)
    res = '{}\t{}\t{}'.format(dict_k[tf], dict_n[tf], pvalue) 
    return res


def gene_binding_count4hypergeometric(in_genes, gene_tf_map):
    '''
    Output:
    a dict of the number genes for binding TFs
    { 'TF1' : 10, 'TF2': 20 }
    '''
    c = Counter()
    input_genes = set(in_genes)
    for i in input_genes:
        if i in gene_tf_map:
            tfs = set(gene_tf_map[i].replace(',', ':').split(':'))
            c.update(tfs)
    return c
##########################################################################################

def obtain_all_TF_biosample_pairs(path_to_all_pairs):
   
    '''iterate all TF_biosample pairs to do TF enrichment for every TF_biosample pair
       during every iteration, permutation is done by the function: "permutation"
    '''
    # ./refs/ENCODE/
    pat = re.compile(r'^(.*)\((.*)\)$')
    path = os.path.join(path_to_all_pairs,'')
    pairs = set()
    tf = None
    sample = None
    for filename in glob.iglob(path + '*.bed'):
        with io.open(filename, 'r', encoding="utf-8") as f:
            line = f.readline()

        tmp = line.split('\t')[-1]
        mat = pat.search(tmp)
        if mat:
            tf = mat.group(1)
            sample = mat.group(2)
        pairs.add((tf,sample,filename,))
    return pairs

##############################################################################
##############################################################################
# -------------------------------------------------------------------------- #
# PCC calculator and network generator based on enriched TFs and input genes #
# -------------------------------------------------------------------------- #


def tf_gene_heatmap(tgt_genes_dict, gene_tf_map, gene_tss_number, outfile):
    '''
    consctruct TF and gene relationship in a heatmap-like format
    actually, in a tab file
    '''
    #output_file = output_file.replace('table','heatmap')

    #folder = './network/'
    #if not os.path.exists(folder):
    #    os.makedirs(folder)

    tfs = {}
    tgt_genes = tgt_genes_dict.keys()
    for gene in tgt_genes:
        try:
            tfs[gene] = gene_tf_map[gene]
        except KeyError:
            pass
    # no genes found have TF peaks
    if len(tfs) == 0:
        print('No TFs found for input genes!') 
        #sys.exit(1)
        return None
    else:
        '''
        #gene => TF1:TF2:TF1,TF2:TF3
        #[(TF1,),(TF2,),(TF1,TF2,),(TF3,)]
        '''
        all_tfs = set()
        gene_tf_relation = dict((k, {})for k in tfs)
        for gene in tfs:
            num = gene_tss_number[gene]
            tf_list = []
            peaks = tfs[gene].split(':')
            for peak in peaks:
                l = peak.split(',')
                l.sort()
                tf_list.append(tuple(l))

            count = Counter()
            for tf in tf_list:
                count.update(tf)
            for i in count:
                count[i] = round(count[i]/float(num),4)
                gene_tf_relation[gene][i] = count[i]
                all_tfs.add(i)
        output = io.open(outfile, 'w', encoding='utf-8')
        ordered_genes = list(gene_tf_relation.keys())
        ordered_genes.sort()
        ordered_genes_symbol = []
        for i in ordered_genes:
            ordered_genes_symbol.append(tgt_genes_dict[i])
        output.write('TF\t{}\n'.format('\t'.join(ordered_genes_symbol)))
        for tf in all_tfs:
            output.write(tf + '\t')
            for gene in ordered_genes:
                if gene != ordered_genes[-1]:
                    try:
                        output.write('{}\t'.format(gene_tf_relation[gene][tf]))
                    except KeyError:
                        output.write('{}\t'.format(0.0))
                else:
                    try:
                        output.write('{}\n'.format(gene_tf_relation[gene][tf]))
                    except KeyError:
                        output.write('{}\n'.format(0.0))
        output.close()    
#############################################################################################
############################################################################################
#####################################################################################################
def merge_one_bed_old(basefile,tf,cell_type,input_bed,fraction=0.67,only_one=True,output_name=None):
    # annotation = "{0}({1})".format(tf, cell_type)
    if not output_name:
        output_name = "{}_{}.bed".format(tf, cell_type).replace(' ', "_").replace('/', '_')
    cmd = 'bedtools coverage -f {} -a {} -b {} > {}.tmp'.format(fraction,basefile,input_bed,output_name)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    output = open(output_name, 'w')
    with open(output_name+'.tmp', 'r') as f:
        for line in f:
            l = line.rstrip('\n').split('\t')
            chrm = l[0]
            start = l[1]
            end = l[2]
            fraction_A = float(l[6])
            # 1   90020   90170   USF1    0  USF1(A549)
            if only_one:
                if fraction_A != 0.0:
                    output.write('{0}\t{1}\t{2}\t{3}\t0\t{3}({4})\n'.format(chrm,start,end,tf,cell_type))
            else:
                output.write('{0}\t{1}\t{2}\t{3}\t0\t{3}({4})\n'.format(chrm,start,end,tf,cell_type))
    output.close()
    removefile(output_name + '.tmp')
    return output_name

def merge_one_bed(basefile,tf,cell_type,input_bed,overlap=100,output_name=None):
    # merge the input bed file to remove duplicates
    cmd = 'bedtools merge -i {0} > {0}.merged'.format(input_bed)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    new_input_bed = input_bed + '.merged'
    if not output_name:
        output_name = "{}_{}.bed".format(tf, cell_type).replace(' ', "_").replace('/', '_')
    cmd = 'bedtools intersect -wo -a {} -b {} > {}.tmp'.format(basefile,new_input_bed,output_name)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    output = open(output_name, 'w')
    with open(output_name+'.tmp', 'r') as f:
        for line in f:
            l = line.rstrip('\n').split('\t')
            chrm = l[0]
            start = l[1]
            end = l[2]
            tf_peak_len = int(l[5]) - int(l[4])
            dns_peak_len = int(end) - int(start)
            overlap_len = int(l[6])
            # 1   90020   90170   USF1    0  USF1(A549)
            if tf_peak_len <= dns_peak_len:
                if tf_peak_len < overlap:
                    if overlap_len >= tf_peak_len:
                        output.write('{0}\t{1}\t{2}\t{3}\t0\t{3}({4})\n'.format(chrm,start,end,tf,cell_type))
                elif tf_peak_len >= overlap:
                    if overlap_len >= overlap:
                        output.write('{0}\t{1}\t{2}\t{3}\t0\t{3}({4})\n'.format(chrm,start,end,tf,cell_type))
            else:
                if overlap_len >= overlap:
                    output.write('{0}\t{1}\t{2}\t{3}\t0\t{3}({4})\n'.format(chrm,start,end,tf,cell_type))
    output.close()
    removefile(output_name + '.tmp')
    removefile(new_input_bed)
    return output_name


def merge_beds(basefile,input_fover,output_file,overlap=100):
    output_dict = {}
    with io.open(input_fover, 'r', encoding="utf-8") as f:
        for line in f:
            l = line.rstrip().split()
            tf = l[0]
            cell_type = l[1]
            input_bed = l[2]
            print((tf,cell_type,input_bed))
            if tf and cell_type and input_bed:
                name =  merge_one_bed(basefile,tf,cell_type,input_bed,overlap=100)
                with io.open(name, 'r', encoding="utf-8") as f:
                    for line in f:
                        l = line.rstrip().split()
                        key = '{0}\t{1}\t{2}'.format(*l)
                        value = l[-1] # TF(Cell_type)
                        if not key in output_dict:
                            output_dict[key] = value
                        else:
                            tmp_set = set(output_dict[key].split(','))
                            tmp_set.add(value)
                            output_dict[key] = ','.join(tmp_set)
                            tmp_set = set()
                removefile(name)
            else:
                print("Error in input fover file!")

    output = open(output_file + '.tmp', 'w')
    for pos in output_dict:
        tf = extract_TFs_from_string(output_dict[pos])
        output.write('{}\t{}\t0\t{}\n'.format(pos,tf,output_dict[pos]))
    output.close()
    cmd = 'bedtools sort -i {0}.tmp > {0}'.format(output_file)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error in sorting!'
    else:
        error_msg = ''
        
    removefile(output_file + '.tmp')

    return output_file

def extract_TFs_from_string(input_string):
    pat = re.compile(r'^(.*)\((.*)\)$')
    tf_set = set()
    for i in input_string.split(','):
        mat = pat.search(i)
        if mat:
            tf = mat.group(1)
            cy = mat.group(2)
            tf_set.add(tf)
    return ','.join(tf_set)

class TFParser(argparse.ArgumentParser):
    '''
    Ex version for argparse , add raise error function .
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()
####################################################################################
#####################################################################################
def parse_args():
    '''
    Read parameters
    '''
    description = "A versatile TF enrichment tool"
    parser = TFParser(description = description)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 4.2')
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### permutation test
    permutate_parser = sub_parsers.add_parser("test", help = "Permutation test for enriched TFs",
        description = "(Permutation test for enriched TFs) Usage: %(prog)s -i input_gene_list -p TF.bed -o output_permutation_results -u upstream(bp) -d downstream(bp) -m empirical")
    permutate_parser.add_argument("-i","--input",action='store',type=str,dest='input_file',
        required=True, help = "input gene list supporting ENSEMBL Gene ID, Entrez Gene ID and Gene Symbol [required]" )
    permutate_parser.add_argument("-m","--method",action='store',type=str,dest='method',
        help = "Method used to calculate P-values (default: %(default)s)",choices=['empirical','kde','hypergeometric'],default='empirical')
    permutate_parser.add_argument("-o","--out-prefix",action="store",type=str,required=True,dest="output_prefix",help="Prefix of output files(s). [required]")
    permutate_parser.add_argument("-s","--out-suffix",action="store",type=str, dest="output_suffix",help="suffix of output files(s). (default: %(default)s)", choices = ['csv','tab'], default='tab',)
    permutate_parser.add_argument("-p","--peak", dest='peak',choices=['ENCODE_merged.bed','H1-hESC.bed','HeLa-S3.bed','HepG2.bed','K562.bed','LCL.bed','GM12878.bed','A549.bed','HEK293.bed','MCF-7.bed', 'ALL'],
            default='ALL', 
            help="input TF ChIP-seq peaks with multiple TFs (default: %(default)s)",)
    permutate_parser.add_argument("--path", dest='path', default='./refs/ENCODE/', 
          help="peak path for all BED files, required if you set --peak as ALL (default: %(default)s)",)
    permutate_parser.add_argument("-u","--upstream", dest='upstream',default=2000,type=str,
        help = "the number of nucleotides upstream of TSS/5' end, (default: %(default)s)",)
    permutate_parser.add_argument("-d","--downstream", dest='downstream',default=2000,type=str,
        help = "the number of nucleotides downstream of TSS/5' end, (default: %(default)s)",)
    permutate_parser.add_argument("-t","--thread", dest='thread',default=1,type=int,
        help = "the number of thread used (default: %(default)s)",)
    permutate_parser.add_argument("-n","--number", dest='number',default=1000,type=int,
        help = "the number of permutations used (default: %(default)s)",)
    permutate_parser.add_argument("--background", dest='background',
        help = "the customized background gene list for permutation test (default: 'Reference gene model')",)
    permutate_parser.add_argument("-b","--biotype", dest='biotype', choices = ['protein_coding','noncoding','pseudogene','all'], default='protein_coding',
        help = "the biotype of input gene list (default: %(default)s)",)
    permutate_parser.add_argument("-l","--log", dest='log_file', action="store", type=str, default='tf_star.json', help = 'record the progress of permutations, {"percent": percentage, "message": "No. permutation(s) processed"}(default: %(default)s)',)
    ######################################################################
    #### merge customized TF ChIP-seq BED files
    merge_parser = sub_parsers.add_parser("merge", help = "Merge input ChIP-seq BED files into one",
        description = "Usage: %(prog)s")
    group = merge_parser.add_mutually_exclusive_group()
    group.add_argument("-i","--input",action='store',type=str,dest='input_file',
      help = "Specify the input one BED file to use, please use --file if you have multiple BED files." )
    group.add_argument("-f","--file",dest='fover',action='store',help = "specify the input file containing the information of BED files, format : TF<space/tab>Cell_type<space/tab>path/to/the BED_file_name")
    merge_parser.add_argument("-b","--base",dest='base_file',default='master.peaks.hg19.bed',help = "specify the base BED file to use (default: %(default)s) [required]")
    merge_parser.add_argument("--tf", dest='tf',
                             help = "specify the TF for the input one BED file.", )
    merge_parser.add_argument("--cell", dest='cell_type',
                             help = "specify the cell type for the input one BED file", )
    merge_parser.add_argument("-o","--output",dest='output',required=True, action='store', 
            help = "specify the output file name", )
    merge_parser.add_argument("-t","--thread", dest='thread',default=1,type=int,
        help = "the number of thread used (default: %(default)s)",)

    args = parser.parse_args()
    if args.sub_command == 'merge': 
        if args.input_file and not args.tf: 
            parser.error('--tf is required when -i is set.')
        if args.input_file and not args.cell_type:
            parser.error('--cell is required when -i is set.')
        else:
            pass
            #parser.print_help()
    if args.sub_command == 'test': 
        if args.peak == 'ALL' and not args.path: 
            parser.error('--peak is required when -p/--peak is set.')
        else:
            pass

    return args

# ------------------------------------
# Main function
# ------------------------------------

def main():
    args = parse_args()
    ref = './refs/' 
    if args.sub_command == 'test':
        input_genes, id_type = check_input_id(args.input_file)
        valid_input_genes_dict = id_converter(input_genes, id_type)
        valid_input_genes = set(valid_input_genes_dict.keys())

        if not args.background:
            gene_map = background_genes_map(id_type, args.biotype)
            background_genes = set(gene_map.keys())
        else:
            background_genes, background_id_type = check_input_id(args.background)
            valid_background_genes_dict = id_converter(background_genes, background_id_type)
            background_genes = set(valid_background_genes_dict.keys())
            if id_type != background_id_type:
                print("gene IDs don't match in your input background genes and target genes!")

        tss_ref_file, gene_tss_number = reference_map(id_type=id_type, biotype=args.biotype, upstream=args.upstream, downstream=args.downstream)

        if args.peak != 'ALL':
            gene_tf_map, overlap_bed = overlap_peaks(tss_ref=tss_ref_file, tf_peak=ref+args.peak)
            # generate Heatmap for NcRG inference (only works for merged datasets)
            # 
            tf_gene_heatmap(valid_input_genes_dict, gene_tf_map, gene_tss_number, args.output_prefix+'.heatmap')
            if args.method == 'hypergeometric':
                log_file = args.log_file
                if log_file:
                    log = io.open(log_file, 'w', encoding='utf-8')

                dict_n = gene_binding_count4hypergeometric(in_genes=background_genes, gene_tf_map=gene_tf_map)
                dict_k = gene_binding_count4hypergeometric(in_genes=valid_input_genes, gene_tf_map=gene_tf_map)
                N = len(valid_input_genes)
                M = len(background_genes)
                # binding_NO_input_genes, binding_NO_total_genes, hypergeometric_pvalue
                results_list = {}
                for tf in dict_k:
                    pvalue = hypergeometric_pvalue(dict_k[tf], M, dict_n[tf], N)
                    results_list[tf] = '{}\t{}\t{}'.format(dict_k[tf], dict_n[tf], pvalue) 
                 
                biosample = os.path.splitext(args.peak)[0] 
                output_result(results_list=results_list, output_prefix=args.output_prefix, output_suffix=args.output_suffix, pairwise=False, biosample=biosample,sampling=False) 
                if log_file:
                    message = 'completed'
                    data = {"percent": int(100), "message": message}
                    log.seek(0)
                    log.write(str(json.dumps(data, ensure_ascii=False)))
                    log.close()

            # empirical sampling (KDE pvalue or empirical pvalue)
            else:
                #print gene_tf_map
                results_list = permutation(background_genes=background_genes, gene_tf_map=gene_tf_map, gene_tss_number=gene_tss_number, tgt_genes=valid_input_genes, n=args.number, threads=args.thread, log_file=args.log_file, method=args.method)
                removefile(tss_ref_file)
                removefile(overlap_bed)

                biosample = os.path.splitext(args.peak)[0] 
                output_result(results_list=results_list, output_prefix=args.output_prefix, output_suffix=args.output_suffix, pairwise=False, biosample=biosample) 
            # tf_gene_heatmap(tgt_genes_dict=input_genes_dict, geneID_tf_map=gene2tf, gene_tss_number=gene_tss_num, output_file=args.output_prefix+'.tab')

        # pair-wise permutation
        else:
            workers = min(args.thread, MAX_WORKERS)
            tf_sample_pairs = obtain_all_TF_biosample_pairs(args.path)
            #print(tf_sample_pairs)
            tf_results = []
            total = len(tf_sample_pairs)
            count = 0

            if args.method == 'hypergeometric':
                log_file = args.log_file
                if log_file:
                    log = io.open(log_file, 'w', encoding='utf-8')

                with ProcessPoolExecutor(max_workers=workers) as executor: 
                    to_do = {}
                    for pair in tf_sample_pairs:
                        _tf, _sample, peak_file = pair
                        #hypergeometric_calculation(tss_ref, tf_peak, background_genes, tgt_genes)
                        future = executor.submit(hypergeometric_calculation, tss_ref_file, peak_file, background_genes, valid_input_genes, _tf)
                        to_do[future] = _tf + ':' + _sample
                    for future in tqdm.tqdm(as_completed(to_do),total=total):
                        count = count + 1
                        if log_file:
                            percent = '{}'.format(int(round(100*count/float(total))))
                            message = '{} permutation(s) processed.'.format(count)
                            data = {"percent": int(percent), "message": message}
                            #log.flush()
                            log.seek(0)
                            # output to JSON
                            log.write(str(json.dumps(data, ensure_ascii=False)))
                        try:
                            result = future.result()
                            tf, biosample = to_do[future].split(':')
                            #print('{}:{}'.format(tf,biosample))
                            #print(other)
                            if result != '':
                                tf_results.append('{}\t{}\t{}'.format(tf,biosample,result))
                        except Exception as exc: 
                            error_msg = 'error: ' + str(exc)
                            print(error_msg)
                if log_file:
                    log.close()
                output_result(results_list=tf_results, output_prefix=args.output_prefix, output_suffix=args.output_suffix, pairwise=True, biosample='', sampling=False) 
            else:
                log_file = args.log_file
                if log_file:
                    log = io.open(log_file, 'w', encoding='utf-8')

                with ProcessPoolExecutor(max_workers=workers) as executor: 
                    to_do = {}
                    for pair in tf_sample_pairs:
                        _tf, _sample, peak_file = pair 
                        # running the pairwise permutation
                        future = executor.submit(pairwise_permutation, tss_ref_file, peak_file, background_genes, gene_tss_number, valid_input_genes, args.number, args.method)
                        to_do[future] = _tf + ':' + _sample
                    for future in tqdm.tqdm(as_completed(to_do),total=total):
                        count = count + 1
                        if log_file:
                            percent = '{}'.format(int(round(100*count/float(total))))
                            message = '{} permutation(s) processed.'.format(count)
                            data = {"percent": int(percent), "message": message}
                            #log.flush()
                            log.seek(0)
                            # output to JSON
                            log.write(str(json.dumps(data, ensure_ascii=False)))
                        try:
                            result = future.result()
                            tf, biosample = to_do[future].split(':')
                            #print('{}:{}'.format(tf,biosample))
                            other = list(result.values())[0]
                            #print(other)
                            tf_results.append('{}\t{}\t{}'.format(tf,biosample,other))
                        except Exception as exc: 
                            error_msg = 'error: ' + str(exc)
                            print(error_msg)

                if log_file:
                    log.close()
                output_result(results_list=tf_results, output_prefix=args.output_prefix, output_suffix=args.output_suffix, pairwise=True, biosample='') 

    elif args.sub_command == 'merge':
        if args.input_file and args.tf and args.cell_type:
            merge_one_bed(basefile=ref+args.base_file,tf=args.tf,cell_type=args.cell_type,input_bed=args.input_file,overlap=100,output_name=args.output)
        elif args.fover:
            merge_beds(basefile=ref+args.base_file,input_fover=args.fover,output_file=args.output,overlap=100)
        else:
            pass


if __name__== '__main__':
    try:
        # tf_star.py test -i genes -p TF.bed -u 2000 -d 2000 -n 1000 -o permutation.result -b protein_coding
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
