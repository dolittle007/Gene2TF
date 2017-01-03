#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'T-Y Wang'
# ------------------------------------
# Python Module
# ------------------------------------
from __future__ import with_statement
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
import numpy as np
from time import sleep
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

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

# ------------------------------------
# Tool functions
# ------------------------------------

def parse_input_file(input_file):
    '''input  : gene list file
       output : a set of gene 
    '''
    input_gene_set = set()
    with open(input_file, 'rU') as f:
        for line in f:
            input_gene_set.add(line.rstrip())
    return input_gene_set

def symbol2gene_map(refgene, biotype):
    '''gene symbols to their corresponding genes from reference gene model  
    '''
    symbol2gene_dict = {}
    # 1   11868   14409   ENST00000456328 gene_name   +   ENSG00000223972 pseudogene
    # 1   30365   30503   NR_036051.1 MIR1302-2   +   100302278   noncoding
    with open(refgene, 'rU') as f:
        for line in f:
            l = line.rstrip().split('\t')
            symbol = l[4]
            gene = l[6]
            bio_type = l[7]
            gene_name = l[4]
            if not gene in symbol2gene_dict:
                if biotype_classify(bio_type, biotype):
                    symbol2gene_dict[gene] = gene_name.upper()
            else:
                pass
    return symbol2gene_dict

def gene2symbol(tgt_genes, refgene):
    '''return a dictionary:
    key: gene ID
    value: symbol
    '''
    gene_to_symbols = symbol2gene_map(refgene, 'all')
    ensembl = re.compile(r'^ENSG\d+$')
    entrez = re.compile(r'^\d+$')
    tgt_genes = set(tgt_genes)
    # checking input IDs
    symbol = False
    output_genes = {}
    for gene in tgt_genes:
        if ensembl.match(gene):
            symbol = False
        elif entrez.match(gene):
            symbol = False
        else:
            symbol = True
    if symbol:
        for key, value in gene_to_symbols.iteritems():
            if value in tgt_genes:
                output_genes[key] = value
        if len(output_genes) == 0:
            print 'Please check you input, only ensembl IDs, entrez IDs and gene symbols are allowed!'
            sys.exit(1)
    else:
        for gene in tgt_genes:
            if gene in gene_to_symbols:
                output_genes[gene] = gene_to_symbols[gene]
    return output_genes

def biotype_classify(src_type, tgt_type):
    biotype = {
            '3prime_overlapping_ncrna' : 'long_noncoding',
            'antisense' : 'long_noncoding',
            'IG_C_gene' : 'protein_coding',
            'IG_C_pseudogene' : 'pseudogene',
            'IG_D_gene' : 'protein_coding',
            'IG_J_gene' : 'protein_coding',
            'IG_J_pseudogene' : 'pseudogene',
            'IG_V_gene' : 'protein_coding',
            'IG_V_pseudogene' : 'pseudogene',
            'lincRNA' : 'long_noncoding',
            'miRNA' : 'short_noncoding',
            'misc_RNA' : 'short_noncoding',
            'Mt_rRNA' : 'short_noncoding',
            'Mt_tRNA' : 'short_noncoding',
            'polymorphic_pseudogene' : 'protein_coding',
            'processed_transcript' : 'long_noncoding',
            'protein_coding' : 'protein_coding',
            'pseudogene' : 'pseudogene',
            'rRNA' : 'short_noncoding',
            'sense_intronic' : 'long_noncoding',
            'sense_overlapping' : 'long_noncoding',
            'snoRNA' : 'short_noncoding',
            'snRNA' : 'short_noncoding',
            'TR_C_gene' : 'protein_coding',
            'TR_D_gene' : 'protein_coding',
            'TR_J_gene' : 'protein_coding',
            'TR_J_pseudogene' : 'pseudogene',
            'TR_V_gene' : 'protein_coding',
            'TR_V_pseudogene' : 'pseudogene',
            'noncoding' : 'noncoding', # NCBI refseq gene
    }
    classification_types = {'protein_coding','pseudogene', 'long_noncoding', 'short_noncoding','noncoding','all'}
    if not tgt_type in classification_types:
        print 'check you input biotype!'
        sys.exit(1)
        return None
    else:
        if biotype[src_type] == tgt_type:
            return True
        elif tgt_type == 'noncoding':
            if biotype[src_type] == 'long_noncoding' or biotype[src_type] == 'short_noncoding':
                return True
        elif tgt_type == 'all':
            return True



def gene2tf_map(refgene, tf_peaks, upstream, downstream, biotype):
    '''construct gene to TFs relationships
    '''
    tmp_folder = './tmp/'
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
    tss_bed_file = tmp_folder + 'tss_upstream_' + str(upstream) + '_downstream_' + str(downstream) + '.bed'

    tss_bed = open(tss_bed_file, 'w')
    stored_region_dict = {}
    gene_tss_number = defaultdict(int)
    with open(refgene, 'rU') as f:
        for line in f:
            l = line.rstrip().split('\t')
            chrm  = l[0]
            start = l[1]
            end   = l[2]
            transcript = l[3]
            strand = l[5]
            gene = l[6]
            bio_type = l[7]
            if strand == '+' and biotype_classify(bio_type, biotype):
                new_start = 0 if int(start) - int(upstream) <=0 else int(start) - int(upstream)
                new_end   = CHRS[chrm] if int(start) + int(downstream) >= CHRS[chrm] else int(start) + int(downstream)
                stored_region = chrm + ':' + strand + ':' + str(new_start) + ':' + str(new_end)
                if not stored_region in stored_region_dict:
                    stored_region_dict[stored_region] = True
                    tss_bed.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrm,new_start,new_end,transcript,'.',strand,gene,bio_type))
                    gene_tss_number[gene] = gene_tss_number[gene] + 1
            elif strand == '-' and biotype_classify(bio_type, biotype):
                new_start = 0 if int(end) - int(downstream) <= 0 else int(end) - int(downstream)
                new_end   = CHRS[chrm] if int(end) + int(upstream) >= CHRS[chrm] else int(end) + int(upstream)
                stored_region = chrm + ':' + strand + ':' + str(new_start) + ':' + str(new_end)
                if not stored_region in stored_region_dict:
                    stored_region_dict[stored_region] = True
                    tss_bed.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrm,new_start,new_end,transcript,'.',strand,gene,bio_type))
                    gene_tss_number[gene] = gene_tss_number[gene] + 1
    tss_bed.close()
    # release the memory
    stored_region_dict = {}
    # intersect TSS bed to TF ChIP-seq bed
    cmd = 'bedtools intersect -wo -a {} -b {} > {}overlap.bed'.format(tss_bed_file,tf_peaks,tmp_folder)
    try: 
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    if error_msg:
        print '*** Error for generating overlapping bed file: \n{}\n{}'.format(err, error_msg)

    # overlapping bed filtering 
    # 1   9873    13873   ENST00000518655 .   +   ENSG00000223972 pseudogene  1   10180   10330   CEBPB,ZBTB33    28  CEBPB(H1-hESC),ZBTB33(HepG2),ZBTB33(K562)   proximal    150
    geneID_tf_map = {}
    with open(tmp_folder + 'overlap.bed', 'rU') as f:
        for line in f:
            l = line.rstrip().split('\t')
            geneID     = l[6]
            tfs        = l[11]
            tss_start  = int(l[1])
            tss_end    = int(l[2])
            peak_start = int(l[9])
            peak_end   = int(l[10])
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
                if not geneID in geneID_tf_map:
                    geneID_tf_map[geneID] = tfs
                else:
                    geneID_tf_map[geneID] = geneID_tf_map[geneID] + ':' + tfs

    removefile(tss_bed_file)
    removefile(tmp_folder + 'overlap.bed') 

    print "gene TF map loaded!"
    # modified Oct. 4th, 2016
    return (geneID_tf_map, gene_tss_number)

def removefile(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

def tf_calc(tgt_genes, geneID_tf_map, gene_tss_number):
    tfs = {}
    for gene in tgt_genes:
        try:
            tfs[gene] = geneID_tf_map[gene]
        except KeyError:
            pass
    # no genes found have TF peaks
    if len(tfs) == 0:
        print 'No TFs found for input genes!' 
        #sys.exit(1)
        return None
    else:
        '''
        #gene => TF1:TF2:TF1,TF2:TF3
        #[(TF1,),(TF2,),(TF1,TF2,),(TF3,)]
        '''
        count_total = Counter()
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
                count[i] = count[i]/float(num)
            count_total.update(count)
        return count_total

def tf_gene_heatmap(tgt_genes_dict, geneID_tf_map, gene_tss_number, output_file):
    '''
    consctruct TF and gene relationship in a heatmap-like format
    actually, in a tab file
    '''
    output_file = output_file.replace('table','heatmap')
    #tmp_folder = './tmp/'
    #if not os.path.exists(tmp_folder):
    #    os.makedirs(tmp_folder)
    tfs = {}
    tgt_genes = set(tgt_genes_dict.keys())
    for gene in tgt_genes:
        try:
            tfs[gene] = geneID_tf_map[gene]
        except KeyError:
            pass
    # no genes found have TF peaks
    if len(tfs) == 0:
        print 'No TFs found for input genes!' 
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
        output = open(output_file, 'w')
        ordered_genes = gene_tf_relation.keys()
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

def aggregate_tf_count(observed, distribution):
    tf_names = observed.keys()
    permutation_count = dict((k, 0) for k in tf_names)
    permutation_p = dict((k, 0) for k in tf_names)
    tf_numbers = dict((k, []) for k in tf_names)
    
    tf_numbers_mean = dict((k, 0) for k in tf_names)
    tf_numbers_std = dict((k, 0) for k in tf_names)
    fold_change = dict((k, 0) for k in tf_names)

    for permutation in distribution:
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
        permutation_p[i] = float(permutation_count[i])/len(distribution)
        tf_numbers_mean[i] = round(np.mean(tf_numbers[i]),4)
        if np.sum(tf_numbers[i]) != 0.0:
            tf_numbers_std[i] = round(np.std(tf_numbers[i]), 4)
            fold_change[i] = round(observed[i]/float(np.mean(tf_numbers[i])), 4)
        else:
            tf_numbers_std[i] = 'NaN'
            fold_change[i] = 'NaN'

        permutation_result[i] = '{}\t{}\t{:.4f}\t{}\t{}\t{}'.format(i, permutation_p[i], observed[i], tf_numbers_mean[i], tf_numbers_std[i], fold_change[i])

    return permutation_result


def permutation(background_genes, gene2tf_map, gene_tss_number, tgt_genes, n, threads, output_prefix, log_file, output_suffix):
    print "Permutation begins!"
    print datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')
    workers = min(MAX_WORKERS, threads)
    #background_genes = gene_map.keys()
    observed_tf_count = tf_calc(tgt_genes, gene2tf_map, gene_tss_number)
    
    tf_count_distribution = []

    count = 0
    if log_file:
        log = io.open(log_file, 'w', encoding='utf-8')

    with ThreadPoolExecutor(max_workers=workers) as executor:
        tf_stats = {}
        for i in xrange(n):
            permut_genes = random.sample(background_genes, len(tgt_genes))
            future = executor.submit(tf_calc, permut_genes, gene2tf_map, gene_tss_number)
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
                log.write(unicode(json.dumps(data, ensure_ascii=False)))
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
                print '*** Error for {}: {}'.format(idx, error_msg)
    if log_file:
        log.close()
    permutation_result = aggregate_tf_count(observed_tf_count, tf_count_distribution)

    print "Permutation accomplished!" 
    print datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')

    output_file = output_prefix + '.' + output_suffix
    with open(output_file, 'w') as f:
        header = 'TF\tpermutation_pvalue\tobserved_count\tpermutation_avg_count\tpermutation_count_std\tfold_change\n'
        if output_suffix == 'csv':
            header = header.replace('\t',',')
        f.write(header)    
        for tf in permutation_result:
            tf_result = permutation_result[tf]
            if output_suffix == 'csv':
                tf_result = tf_result.replace('\t',',')
            f.write('{}\n'.format(tf_result))

    return output_file

def merge_one_bed(basefile,tf,cell_type,input_bed,fraction=0.67,only_one=True,output_name=None):
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
    with open(output_name+'.tmp', 'rU') as f:
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

def merge_beds(basefile,input_fover,output_file,fraction=0.67):
    output_dict = {}
    with open(input_fover, 'rU') as f:
        for line in f:
            l = line.rstrip().split()
            tf = l[0]
            cell_type = l[1]
            input_bed = l[2]
            if tf and cell_type and input_bed:
                name = merge_one_bed(basefile,tf,cell_type,input_bed,fraction=fraction,only_one=True)
                with open(name, 'rU') as f:
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
                print "Error in input fover file!"

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

def parse_args():
    '''
    Read parameters
    '''
    description = "TF enrichment tool"
    parser = TFParser(description = description)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2')
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### permutation test
    permutate_parser = sub_parsers.add_parser("test", help = "Permutation test for enriched TFs",
        description = "(Permutation test for enriched TFs) Usage: %(prog)s -i input_gene_list -r ref_gene_model -p TF.encode.bed -o output_permutation_results -u upstream(bp) -d downstream(bp)")
    permutate_parser.add_argument("-i","--input",action='store',type=str,dest='input_file',
        required=True, help = "input gene list supporting ENSEMBL Gene ID, Entrez Gene ID and Gene Symbol [required]" )
    permutate_parser.add_argument("-r","--refgene",action='store',type=str,dest="ref_gene_model",
        required=True,help='Reference gene model in bed format. [required]',)
    permutate_parser.add_argument("-o","--out-prefix",action="store",type=str,required=True,dest="output_prefix",help="Prefix of output files(s). [required]")
    permutate_parser.add_argument("-s","--out-suffix",action="store",type=str, dest="output_suffix",help="suffix of output files(s). (default: %(default)s)", choices = ['csv','tab'], default='tab',)
    permutate_parser.add_argument("-p","--peak", dest='peak',choices=['encode_tf_phase3.bed','H1-hESC.encode.bed','HeLa-S3.encode.bed','HepG2.encode.bed','K562.encode.bed','LCL.encode.bed','GM12878.encode.bed','A549.encode.bed','HEK293.encode.bed','MCF-7.encode.bed'],
            default='encode_tf_phase3.bed', 
            help="input merged TF ChIP-seq peaks with multiple TFs (default: %(default)s)",)
    permutate_parser.add_argument("-u","--upstream", dest='upstream',default=2000,type=int,
        help = "the number of nucleotides upstream of TSS/5' end, (default: %(default)s)",)
    permutate_parser.add_argument("-d","--downstream", dest='downstream',default=2000,type=int,
        help = "the number of nucleotides downstream of TSS/5' end, (default: %(default)s)",)
    permutate_parser.add_argument("-t","--thread", dest='thread',default=1,type=int,
        help = "the number of thread used (default: %(default)s)",)
    permutate_parser.add_argument("-n","--number", dest='number',default=1000,type=int,
        help = "the number of permutations used (default: %(default)s)",)
    permutate_parser.add_argument("--background", dest='background',
        help = "the customized background gene list for permutation test (default: 'Reference gene model')",)
    permutate_parser.add_argument("-b","--biotype", dest='biotype', choices = ['protein_coding','pseudogene','long_noncoding','short_noncoding','noncoding','all'], default='protein_coding',
        help = "the biotype of input gene list (default: %(default)s), 'long_noncoding','short_noncoding' only work for ENSEMBL gene model",)
    permutate_parser.add_argument("-l","--log", dest='log_file', action="store", type=str, default='tf_star.json', help = 'record the progress of permutations, {"percent": percentage, "message": "No. permutation(s) processed"}(default: %(default)s)',)

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
    return args

# ------------------------------------
# Main function
# ------------------------------------

def main():
    args = parse_args()
    ref = './refs/' 
    if args.sub_command == 'test':
        # permutation test
        gene2tf, gene_tss_num = gene2tf_map(refgene=ref+args.ref_gene_model, tf_peaks=ref+args.peak, upstream=args.upstream, downstream=args.downstream, biotype=args.biotype)
        input_genes_dict = gene2symbol(parse_input_file(args.input_file),refgene=ref+args.ref_gene_model)
        input_genes = set(input_genes_dict.keys())

        tf_gene_heatmap(tgt_genes_dict=input_genes_dict, geneID_tf_map=gene2tf, gene_tss_number=gene_tss_num, output_file=args.output_prefix+'.tab')
        if not args.background:
            gene_map = symbol2gene_map(refgene=ref+args.ref_gene_model,biotype=args.biotype)
            background_genes = set(gene_map.keys())
        else:
            background_genes_dict = gene2symbol(parse_input_file(args.background),refgene=ref+args.ref_gene_model)
            background_genes = set(background_genes_dict.keys())

        permutation(background_genes=background_genes, gene2tf_map=gene2tf, gene_tss_number=gene_tss_num, tgt_genes=input_genes, n=args.number, threads=args.thread, output_prefix=args.output_prefix, log_file=args.log_file, output_suffix=args.output_suffix)

    elif args.sub_command == 'merge':
        if args.input_file and args.tf and args.cell_type:
            merge_one_bed(basefile=ref+args.base_file,tf=args.tf,cell_type=args.cell_type,input_bed=args.input_file,fraction=0.67,only_one=True,output_name=args.output)
        elif args.fover:
            merge_beds(basefile=ref+args.base_file,input_fover=args.fover,output_file=args.output,fraction=0.67)
        else:
            pass


if __name__== '__main__':
    try:
        # tf_star.py test -i genes -r Homo_sapiens.GRCh37.ensembl.bed -p TF.encode.bed -u 2000 -d 2000 -n 1000 -o permutation.result -b protein_coding
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
