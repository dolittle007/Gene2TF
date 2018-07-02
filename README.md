
Introduction
------------
Gene2TF: a versatile tool for detecting transcription factors underlying genes or variants in diseases

Prerequisites
------------
The software is developed and tested in Linux and Max OS environments.

You need [BEDTools](https://github.com/arq5x/bedtools2),  [Python 3.5](https://www.python.org/), and some python packages 

* [numpy](http://www.numpy.org/)
* [tqdm](https://pypi.python.org/pypi/tqdm)
* [argparse](https://pypi.python.org/pypi/argparse)

to run Gene2TF.

Getting Soure Code
------------------
	git clone git://github.com/dolittle007/gene2tf.git
	cd gene2tf
Running Gene2TF 

Subcommand
-----------------
  test         Permutation test for enriched TFs
  merge        Merge input ChIP-seq BED files into one

## Permutation test for enriched TFs

gene2tf.py test [-h] -i INPUT_FILE [-m {empirical,kde,hypergeometric}]
                       -o OUTPUT_PREFIX [-s {csv,tab}]
                       [-p {ENCODE_merged.bed,H1-hESC.bed,HeLa-S3.bed,HepG2.bed,K562.bed,LCL.bed,GM12878.bed,A549.bed,HEK293.bed,MCF-7.bed,ALL}]
                       [--path PATH] [-u UPSTREAM] [-d DOWNSTREAM] [-t THREAD]
                       [-n NUMBER] [--background BACKGROUND]
                       [-b {protein_coding,noncoding,pseudogene,all}]
                       [-l LOG_FILE]


Input data
-----------------
Gene2TF will calculate the TF enrichment results from input genes.
Some support data is needed, that needs to be set up prior Gene2TF execution.

The gist of Gene2TF input is:
- database (an example is 
- 
-  results (such aswhich are just randomly generated)
