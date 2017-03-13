# TF-STAR

TF-STAR: a versatile tool for identifying key transcription factors for human diseases through analysis of the regulatory regions of disease genes or associated variants 



## Introduction



## Prerequisites

The software is developed and tested in Linux and Max OS environments.

You need  [Python 2.7](https://www.python.org/), and some python packages 

* [numpy](http://www.numpy.org/)
* [concurrent.futures](https://pypi.python.org/pypi/futures)
* [tqdm](https://pypi.python.org/pypi/tqdm)
* [argparse](https://pypi.python.org/pypi/argparse)

to run TF-STAR.

## Project Layout


## Input data
TF-STAR will calculate the TF enrichment results from input genes.
Some support data is needed, that needs to be set up prior TF-STAR execution.

The gist of TF-STAR input is:
- A Transcriptome Prediction Model database (an example is [here](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/DGN-WB_0.5.db))
- A file with the covariance matrices of the SNPs within each gene model (such as [this one](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/covariance.DGN-WB_0.5.txt.gz))
- GWAS results (such as [these](https://s3.amazonaws.com/imlab-open/Data/MetaXcan/sample_data/GWAS.tar.gz), which are just randomly generated)
