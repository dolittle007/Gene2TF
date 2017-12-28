#!/usr/bin/env python
#-*- coding: utf-8 -*-

##network analysis functions for TF-STAR
from __future__ import division, with_statement
import pandas as pd
import numpy as np
import re
import os
import sys
import random
import time
from collections import defaultdict
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import io
from pypathway import Reactome, GO, KEGG, ORA
from pypathway import IdMapping, GMTUtils
from helper import check_input_id, background_genes_map, id_converter, removefile
from helper import correct_pvalues_for_multiple_testing 

##############################################################################
# -------------------------------------------------------------------------- #
# Network evaluation (the whole network or modules in the network            #
# -------------------------------------------------------------------------- #
def read_net(infile,delimiter='\t',weight=True):
    '''input a weighted network to a networkx object
    '''
    if weight:
        g = nx.read_weighted_edgelist(infile, delimiter=delimiter)
    else:
        g = nx.read_edgelist(infile, delimiter=delimiter)
    return g

def edges_in_ref(tgt, ref):
    '''tgt: target network
       ref: reference network
    '''
    count = 0
    for i in tgt.edges:
        if ref.has_edge(*i):
            count = count + 1
    return count

def net_eval(src_net, ref='./refs/InBio_Map_core.tab', number=1000):
    '''input: src network,
              ref network
              permutation number
       output: empirical p-value based on randomized networks
    '''
    # modify for eval_modules.py 
    #src_net = read_net(src)
    ref_net = read_net(ref)

    if (nx.number_of_nodes(src_net)<=2):
        two_nodes_edge = edges_in_ref(src_net, ref_net)
        if two_nodes_edge == 1:
            return 0.0
        else:
            return 1.0

    src_edge_no = edges_in_ref(src_net, ref_net)
    count = 0 
    for j in range(number):
        src = list(src_net.nodes())
        perms = int(len(src)/2) + 1
        for i in range(1, perms+1):
            random.shuffle(src)
        mapping = dict(zip(src_net.nodes(), src))
        random_net = nx.relabel_nodes(src_net, mapping)

        random_edge_no = edges_in_ref(random_net, ref_net)
        if random_edge_no >= src_edge_no:
            count = count + 1
    if not count:
        return '<{}'.format(1.0/number)
    else:
        return float(count)/number 

#############################################################################
#---------------------------------------------------------------------------#
# Run a Gene Ontology Enrichment Analysis (GOEA) for modules in the network #
#---------------------------------------------------------------------------#
#module_name  gene1,gene2,gene3  P-value

def run_ora(in_genes, background_genes, id_type, database, fig='enrichment.png', plot=True):
    '''in_genes: a list of genes
       id_type: ensembl, entrez, symbol
       database: kegg, reactome, go.bp, go.cc, go.mf, go.all
    '''
    # id_type is 'ensembl', using symbol instead via IdMapping 
    if id_type == 'ensembl':
        id_maps_in = IdMapping.convert(input_id=in_genes, organism='hsa', source='ENSEMBL', target='SYMBOL')
        id_maps_back = IdMapping.convert(input_id=background_genes, organism='hsa', source='ENSEMBL', target='SYMBOL')
        in_genes = [x[1][0] for x in id_maps_in if x[1]]
        background_genes = [x[1][0] for x in id_maps_back if x[1]]
        id_type = 'symbol'

    # reactome.symbol.gmt; reactome.entrez.gmt
    gmt_file = '{}.{}.gmt'.format(database, id_type)
    gmt = GMTUtils.parse_gmt_file('./refs/gmt/' + gmt_file)

    res = ORA.run(in_genes, background_genes, gmt)
    df = res.table
    df.index = df['name'] 

    df_f = df[(df['NDE'] >= 1)]
    df_f.drop(['fdr'], axis=1, inplace=True)
    # calculate FDR based on all the p-values
    df_f.loc[:,'fdr'] = correct_pvalues_for_multiple_testing(df_f['p-value'])  

    #df_filtered = df[(df['NDE'] >= 1) & (df['p-value'] < 0.05)]
    df_filtered = df_f

    
    # if no enrichment found, it is unnecessary to plot 
    if df_filtered.empty:
        return df_filtered
    if plot:
        # for bar plot only 
        df_filtered2 = df_f[(df_f['NDE'] >= 2) & (df_f['fdr'] < 0.05)]

        df_filtered2.loc[:,'log_fdr'] = -np.log10(df_f['fdr'])
        df_draw = df_filtered2.sort_values("log_fdr", ascending=False)
        df_draw = pd.DataFrame(df_draw['log_fdr'])
        # plot only the top ten terms
        df_fig = df_draw.iloc[0:10,]

        #fig = plt.figure(figsize=(4, 5), dpi=100)
        ax = df_fig.plot(kind='barh', figsize=(18,10), legend=False, fontsize=13)
        ax.set(xlabel="$-Log_{10}(FDR)$")
        ax.xaxis.label.set_size(14)
        ax.set(ylabel="")
        plt.title("Enrichment result", size=20)
        try:
            plt.tight_layout()
        except ValueError:
            pass
        finally:
            plt.savefig(fig, dpi=300)
            plt.close()

    # 'name', 'mapped', 'NDE', 'p-value', 'fdr', 'log_fdr', 'DE'
    return df_filtered

def output_enrichments_info_to_excel(df_map,xls_name):
    '''input: multiple ORA results {'module name': df}
       output: one Excel
    '''
    count = 0
    columns = ['name','mapped','NDE', 'p-value', 'fdr', 'DE']
    writer = pd.ExcelWriter(xls_name)
    for name in df_map:
        if not df_map[name].empty:
            count = count + 1
            df_map[name].to_excel(writer, name, columns=columns, index=False)
    writer.save()

    if count == 0:
        return None
    else:
        return xls_name

def enrich_modules(infile, database, bk_genes, id_type):
    '''enrichment analysis for modules
    '''
    folder = './{}'.format(database)
    if not os.path.exists(folder):
        os.makedirs(folder)
    df_map = {}
    with io.open(infile, 'r') as f:
        for line in f:
            name, genes, pvalue = line.rstrip('\n').split()
            valid_input_genes_dict = id_converter(set(genes.split(',')), id_type)
            valid_input_genes = set(valid_input_genes_dict.values())
            df = run_ora(in_genes=valid_input_genes, background_genes=bk_genes, id_type=id_type, database=database,fig='{}/{}.png'.format(folder, name))
            df_map[name] = df

    output_enrichments_info_to_excel(df_map, xls_name='{}/{}.xlsx'.format(folder, database))


##############################################################################
# -------------------------------------------------------------------------- #
# PCC calculator and network generator based on enriched TFs and input genes #
# -------------------------------------------------------------------------- #


def tf_gene_heatmap(tgt_genes, gene_tf_map, gene_tss_number, output_file):
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
    
if __name__ == '__main__':
    #input_ref_net()
    #print(net_eval('./input.net', number=10))
    input_genes, id_type = check_input_id('./B_cell.up.txt')
    valid_input_genes_dict = id_converter(input_genes, id_type)
    if id_type == 'symbol':
        valid_input_genes = set(valid_input_genes_dict.values())
    else:
        valid_input_genes = set(valid_input_genes_dict.keys())

    gene_map = background_genes_map(id_type, 'protein_coding')
    if id_type == 'symbol':
        background_genes = set(gene_map.values())
    else:
        background_genes = set(gene_map.keys())

    df1 = run_ora(in_genes=valid_input_genes, background_genes=background_genes, id_type=id_type, database='reactome',fig='reactome.png')
    df2 = run_ora(in_genes=valid_input_genes, background_genes=background_genes, id_type=id_type, database='go.bp')

    df_map = {'reactome':df1, 'go.bp':df2}

    output_enrichments_info_to_excel(df_map,xls_name='B_cell.enrichment.xlsx')
