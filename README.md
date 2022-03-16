# pykallisto
A Python API for RNA-seq alignment, quantification normalization via Kallisto. PyKallisto can take a list of raw fastq files in a compressed or uncompressed format, and generate a cell x gene expression matrix and/or
 a cell x transcript equivalent class count matrix, quantify data from many common technologies such as 10x, indrops and dropseq, process feature barcoding data, and obtain QC reports for single-cell rna-seq data. Under the hood, PyKallisto is entirely reliant on `kallisto`, `bustools`, and `kallisto | bustools` from the [Patcher Lab](https://pachterlab.github.io). 

 All thanks for the technologies and beautiful command line tools goes to them. All I did was wrap them in a Python API for easy handling of data in larger pipelines, to allow for pure Python single-cell analysis in a single line of code.

## Usage
`pykallisto` has two classes, `Kallisto` and `KallistoBus`. The former wraps `kallisto`, and the latter wraps the `kb | tools` from the Patcher Lab. 

For both tools, the arguments list are the same as the extended argument list in the command line version. That is, instead of `-o/--output-dir` for the output directory, we pass `output_dir=/path/to/dir` in the Python classes. 

### `pykallisto.Kallisto` 
Example:
```python
from pykallisto import Kallisto 

data = Kallisto(files=['fastq1.fasta', 'fastq2.fasta', 'fastq3.fasta'])

# Generate the index file for pseudoalignment 
data.index(index='index', output_dir='index_directory/')

# Quantify transcripts in parallel on 4 threads
data.quant(output_dir='results/', threads=4)
```

Methods:

    Computes equivalence classes for reads and quantifies abundances (deprecated).

    Arguments:

    index: Filename for the kallisto index to be used for pseudoalignment. If Kallisto class was initialized with index file name, this is optional
    files: List of fastq file paths, or path to folder containing fastq files
    output_dir: Directory to write output to
    umi: First file in pair is a UMI file
    batch: Process files listed in FILE
    quant: Quantify using EM algorithm (only in batch mode)
    bus: Output a BUS file
    single: Quantify single-end reads
    fragment_length: Estimated average fragment length
    sd: Estimated standard deviation of fragment length (default: -l, -s values are estimated from paired end data, but are required when using --single unless outputting a BUS file via --bus)
    fr_stranded: Strand specific reads, first read forward
    rf_stranded: Strand specific reads, first read reverse
    num: Output number of read in BUS file flag column (only with --bus)
    threads: Number of threads to use (default: 1)Parameters:

