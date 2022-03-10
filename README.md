# pykallisto
A Python API for RNA-seq alignment, quantification normalization via Kallisto. PyKallisto can take a list of raw fastq files in a compressed or uncompressed format, and generate a cell x gene expression matrix and/or
 a cell x transcript equivalent class count matrix, quantify data from many common technologies such as 10x, indrops and dropseq, process feature barcoding data, and obtain QC reports for single-cell rna-seq data. Under the hood, PyKallisto is entirely reliant on `kallisto`, `bustools`, and `kallisto | bustools` from the [Patcher Lab](https://pachterlab.github.io). 

 All thanks for the technologies and beautiful command line tools goes to them. All I did was wrap them in a Python API for easy handling of data in larger pipelines, to allow for pure Python single-cell analysis in a single line of code.

## Usage
Initialize a Kallisto class object via 
```python
from pykallisto import Kallisto 

data = Kallisto(files=['fastq1.fasta', 'fastq2.fasta', 'fastq3.fasta'])

# Generate the index file for pseudoalignment 
data.index(index='index', output_dir='index_directory/')

# Quantify transcripts in parallel on 4 threads
data.quant(output_dir='results/', threads=4)
```
