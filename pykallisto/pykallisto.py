from multiprocessing.sharedctypes import Value
import os 
import pathlib 
from typing import *
import subprocess

class Kallisto:
    def __init__(
        self, 
        files: Union[List[str], str]=None,
        index_file=None 
    ):  
        self.files = self._validate_files(files)
        self._index_file = (index_file if index_file else None) 

    def _validate_files(files):
        """
        Handles files being a list of fastq files, or the path to a folder containing fastq files.

        Input:
        files: Union[str, List[str]]: List of files or path to folder 

        Returns:
        List[str]: List containing absolute paths to all fastq files
        """

        if isinstance(files, str):
            files = [os.path.join(files, f) for f in os.listdir(files)]
        elif isinstance(files, list):
            files = [f.strip() for f in files]
        else:
            raise ValueError("files must be path to files or list of files")

        return ' '.join(files)

    def _validate_input(self, files, index):
        """
        Validates files input and index input, since the Kallisto class can be initialized with either or have either passed for maximum flexibility

        Input:
        files: Union[str, List[str]]: List of files or path to folder 
        index: Index file name 

        Returns:
        List[str], str: List containing fastq file paths, and index file name 
        """
        # If no files are passed and class wasn't initialized with any 
        if not files and not self.files:
            raise ValueError("One of self.files or files argument must not be None.")

        if files and self.files:
            print('Warning: Class initialized with files list & files list passed. Using passed list.')

        if self._index_file and index:
            print('Warning: Index file generated and index file argument passed, using passed argument.')

        if not self._index_file and not index:
            raise ValueError("Index files not generated. Run self.index() or pass the path to the index file.")

        files = self._validate_files(files)

        return files, index 

    def index(
        self,
        index: str=None,
        files: Union[str, List[str]]=None,
        kmer_size: int=None,
        make_unique: bool=False 
    ):
        """
        Builds kallisto index via pseudoalignment

        Arguments:
        index: Path to index file. If Kallisto class was initialized with an index filename, this argument is optional.
        files: List of fastq files or path to folder
        kmer_size: k-mer (odd) length (default: 31, max value: 31)
        make_unique: Replace repeated target names with unique names

        Returns:
        None
        """
        files, index = self._validate_input(index, files)

        command = (
            "kallisto index "
            f"--index={index} "
            f"{f'--kmer-size={kmer_size} ' if kmer_size else ''}"
            f"{'--make-unique ' if make_unique else ''}"
            f"{files}"
        )

        os.system(
            command 
        )

    def quant(
        self, 
        output_dir: str, 
        index: str=None,
        files: Union[List[str], str]=None,
        bias: bool=False,
        bootstrap_samples: int=None,
        seed: int=42, 
        plaintext: bool=False,
        fusion: bool=False,
        single: bool=False,
        single_overhang: bool=False,
        fr_stranded: bool=False,
        rf_stranded: bool=False,
        sd: bool=False,
        threads: int=1, 
        pseudobam: bool=False,
        genomebam: bool=False,
        gtf: str = None,
        chromosomes: str=None,
        verbose=False
    ):
        """
        Computes equivalence classes for reads and quantifies abundances.

        Arguments:
        output_dir: Directory to write output to
        index: Filename for kallisto index to be used. If Kallisto class was initialized with a index filename, this argument is optional
        files: List of fastq file paths, or path to folder containing fastq files 
        bias: Perform sequence based bias correction
        bootstrap_samples: Number of bootstrap samples (default: 0)
        seed: Seed for the bootstrap sampling (default: 42)
        plaintext: Output plaintext instead of HDF5
        fusion: Search for fusions for Pizzly
        single: Quantify single-end reads
        single-overhang: Include reads where unobserved rest of fragment is predicted to lie outside a transcript
        fr_stranded: Strand specific reads, first read forward
        rf_stranded: Strand specific reads, first read reverse
        fragment_length: Estimated average fragment length
        sd: Estimated standard deviation of fragment length (default: -l, -s values are estimated from paired end data, but are required when using --single)
        threads: Number of threads to use (default: 1)
        pseudobam: Save pseudoalignments to transcriptome to BAM file
        genomebam: Project pseudoalignments to genome sorted BAM file
        gtf: GTF file for transcriptome information (required for --genomebam)
        chromosomes: Tab separated file with chromosome names and lengths(optional for --genomebam, but recommended)
        verbose: Print out progress information every 1M proccessed reads

        Returns:
        None
        """

        files, index = self._validate_input(index, files)

        if not files:
            files = self.files
        index = (index if index else self._index_file)

        files = ' '.join(files)

        command = (
            f"kallisto quant " 
            f"--index={index} "
            f"--output={output_dir} "
            f"{'--bias ' if bias else ''}"
            f"{f'--bootstrap_samples={bootstrap_samples} ' if bootstrap_samples else ''}"
            f"--seed={seed} "
            f"{'--plaintext ' if plaintext else ''}"
            f"{'--fusion ' if fusion else ''}"
            f"{'--single ' if single else ''}"
            f"{'--single_overhang ' if single_overhang else ''}"
            f"{'--fr_stranded ' if fr_stranded else ''}"
            f"{'--rf_stranded ' if rf_stranded else ''}"
            f"{'--sd ' if sd else ''}"
            f"{'--threads ' if threads else ''}"
            f"{'--pseudobam ' if pseudobam else ''}"
            f"{'--genomebam ' if genomebam else ''}"
            f"{'--gtf ' if gtf else ''}"
            f"{'--chromosomes ' if chromosomes else ''}"
            f"{'--verbose ' if verbose else ''}"
            f"{files}"
        )

        os.system(
            command 
        )

    def bus(
        self,
        output_dir: str,
        files: Union[List[str], str]=None,
        technology: str=None,
        index: str=None,
        list: bool=False,
        batch: str=None,
        threads: int=1,
        bam: bool=False,
        num: bool=False,
        tag: str=None,
        fr_stranded: bool=False,
        rf_stranded: bool=False,
        unstranded: bool=False,
        paired: bool=False,
        genomebam: bool=False,
        chromosomes: bool=False,
        gtf: bool=False,
        verbose: bool=False,
    ) -> None:


        """
        Generates BUS files for single-cell sequencing
        
        output_dir: Directory to write output to
        files: List of fastq file paths or path to folder containing fastq files 
        technology: Single-cell technology used 
        index: Filename for the kallisto index to be used for pseudoalignment. If Kallisto class was initialized with index file name, this is optional
        list: List all single-cell technologies supported
        batch: Process files listed in FILE
        threads: Number of threads to use (default: 1)
        bam: Input file is a BAM file
        num: Output number of read in flag column (incompatible with --bam)
        tag: sequence to identify UMI reads for certain technologies
        fr_stranded: Strand specific reads for UMI-tagged reads, first read forward
        rf_stranded: Strand specific reads for UMI-tagged reads, first read reverse
        unstranded: Treat all read as non-strand-specific
        paired: Treat reads as paired
        genomebam: Project pseudoalignments to genome sorted BAM file
        gtf: GTF file for transcriptome information (required for --genomebam)
        chromosomes: Tab separated file with chromosome names and lengths (optional for --genomebam, but recommended)
        verbose: Print out progress information every 1M proccessed reads

        Returns:
        None 
        """

        files, index = self._validate_input(files, index)

        files = ' '.join(files)
        command = (
            f"kallisto bus "
            f"--index={index} "
            f"--output-dir={output_dir} "
            f"--technology={technology} "
            f"{'--list ' if list else ''}"
            f"{f'--batch={batch} ' if batch else ''}"
            f"--threads={threads} "
            f"{'--bam ' if bam else ''}"
            f"{'--num ' if num else ''}"
            f"{f'--tag={tag} ' if tag else ''}"
            f"{'--fr-stranded ' if fr_stranded else ''}"
            f"{'--rf-stranded ' if rf_stranded else ''}"
            f"{'--unstranded ' if unstranded else ''}"
            f"{'--paired ' if paired else ''}"
            f"{'--genomebam ' if genomebam else ''}"
            f"{'--chromosomes ' if chromosomes else ''}"
            f"{'--gtf ' if gtf else ''}"
            f"{'--verbose ' if verbose else ''}"
            f"{files}"
        )

        os.system(
            command 
        )

    def pseudo(
        self, 
        output_dir: str,
        files: Union[List[str], str]=None,
        index: str=None,
        umi: bool=False,
        batch: str=None,
        single: bool=False,
        fragment_length: float=None,
        sd: float=None,
        threads: int=1,
    ) -> None:
        """
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
        threads: Number of threads to use (default: 1)

        Returns:
        None
        """

        files, index = self._validate_input(index, files)

        command = ("kallisto pseudo "
            f"--index={self._index_file} ",
            f"--output-dir={output_dir} ",
            f"{'--umi ' if umi else ''}",
            f"{f'--batch={batch} ' if batch else ''}"
            f"{'--single ' if single else ''}"
            f"{f'--fragment-length={fragment_length} ' if fragment_length else ''}"
            f"{f'--sd={sd} ' if sd else ''}"
            f"--threads={threads} "
            f"{files}"
        )

        os.system(
            command 
        )

    def h5dump(
        self, 
        output_dir: str
    ):

        command = (
            "kallisto h5dump "
            f"--output-dir={output_dir}"
        )

        os.system(
            command
        )

    def merge(
        self,
        index: str,
        output_dir: str 
    ):
        command = (
            "kallisto merge "
            f"--index={index} "
            f"--output-dir={output_dir}"
        )

        os.system(
            command 
        )

    def inspect(
        self,
        gfa: str=None,
        gtf: str=None,
        bed: str=None,
    ):
        command = (
            "kallisto inspect "
            f"{f'--gfa={gfa} ' if gfa else ''}"
            f"{f'--gtf={gtf} ' if gtf else ''}"
            f"{f'--bed={bed} ' if bed else ''}"
        )

        os.system(
            command
        )

class KallistoBus:
    def __init__(
        self,
        files: List[str],
    ) -> None:
        self.files = files 
        
    def info(self):
        os.system(
            "kb info"
        )

    def ref(
        self, 
        index: str,
        t2g: str, 
        fasta: str=None,
        gtf: str=None,
        feature: str=None,
        tmp: str=None,
        keep_temp: bool=False,
        verbose: bool=False,
        include_attribute: Dict[str, str]=None,
        exclude_attribute: Dict[str, str]=None,
        d: str=None,
        k: str=None,
        workflow: str='standard',
        overwrite: bool=False,
        kallisto: str=None,
        bustools: str=None,
        f2: str=None,
        c1: str=None,
        c2: str=None
    ) -> None:

        if fasta and not d:
            raise ValueError("Error, -f1 FASTA argument is only valid with -d passed")

        if d and not fasta:
            raise ValueError("Error, -f1 FASTA argument required with optional argument -d passed.")

        if workflow == 'lamano' or workflow == 'nucleus':
            if not f2 or not c1 or not c2:
                raise ValueError("Error, optional arguments -f2, -c1, -c2 are required with lamano or nucleus workflows.")

        if include_attribute:
            include_attribute_str = [f'--include-attributes={key}:{val} ' for key, val in include_attribute]
            include_attribute_str = ' '.join(include_attribute)

        if exclude_attribute:
            exclude_attribute_str = [f'--exclude-attributes={key}:{val} ' for key, val in exclude_attribute]
            exclude_attribute_str = ' '.join(exclude_attribute_str)
        
        if feature and not workflow == 'kite':
            raise ValueError("Error: Positional argument 'feature' only valid with 'kite' workflow. ")
        
        fasta = ' '.join(fasta)

        command = (
            "kb ref "
            f"-i {index} "
            f"-g {t2g} "
            f"{f'-f1={fasta} ' if fasta else ''}"
            f"{'--tmp ' if tmp else ''}"
            f"{'--keep-tmp ' if keep_temp else ''}"
            f"{'--verbose ' if verbose else ''}"
            f"{include_attribute_str if include_attribute else ''}"
            f"{exclude_attribute_str if exclude_attribute else ''}"
            f"{f'-d {d} ' if d else ''}"
            f"{f'-k {k} ' if k else ''}"
            f"{f'workflow={workflow} ' if workflow else ''}"
            f"{'--overwrite ' if overwrite else ''}"
            f"{f'--kallisto={kallisto} ' if kallisto else ''}"
            f"{f'--bustools={bustools} ' if bustools else ''}"
            f"{fasta} {gtf} {feature if feature else ''}"
        )

        os.system(
            command 
        ) 

    def count(
        self,
        index: str,
        t2g: str,
        technology: str,
        fastqs: List[str]=None,
        tmp: bool=False,
        keep_tmp: bool=False,
        verbose: bool=False,
        out: str=None,
        whitelist: str=None,
        threads: int=1,
        memory: str=None,
        strand: str=None,
        workflow: str=None,
        em: bool=False,
        umi_gene: bool=False,
        mm: bool=False,
        tcc: bool=False,
        filter: str=None,
        filter_threshold: int=None,
        overwrite: bool=False,
        dry_run: bool=False,
        loom: bool=False,
        h5ad: bool=False,
        cellranger: bool=False,
        gene_names: bool=False,
        report: bool=False,
        kallisto: str=None,
        bustools: str=None,
        c1: str=None,
        c2: str=None,
        parity: str=None,
        fragment_l: int=None,
        fragment_s: int=None,
    ) -> None:
        if workflow == 'lamanno' or workflow == 'nucleus' and not c1 and not c2:
            raise ValueError("Optional flags -c1 and -c2 are required when using lamanno or nucleus workflows.")
        
        if parity or fragment_l or fragment_s and (not technology == 'BULK' or not technology == 'SMARTSEQ2'):
            raise ValueError("--parity, --fragment-s and --fragment-l optional parameters are only valid for technologies: BULK or SMARTSEQ2")
        
        if not fastqs:
            fastqs = self.files
        
        if fastqs and self.files:
            print('Warning: KallistoBus object initialized with fastq file list and files list passed. Using arguments past list.')

        files = ' '.join(files)

        command = (
            "kb count "
            f"-i {index} "
            f"-g {t2g} "
            f"-x {technology} "
            f"{'--tmp ' if tmp else ''}"
            f"{'--keep-tmp ' if keep_tmp else ''}"
            f"{'--verbose ' if verbose else ''}"
            f"{f'-o {out} ' if out else ''}"
            f"{f'-w {whitelist} ' if whitelist else ''}"
            f"{f'-t {threads} ' if threads else ''}"
            f"{f'-m {memory} ' if memory else ''}"
            f"{f'--strand={strand} ' if strand else ''}"
            f"{f'--workflow={workflow} ' if workflow else ''}"
            f"{'--em ' if em else ''}"
            f"{'--umi-gene ' if umi_gene else ''}"
            f"{'-mm ' if mm else ''}"
            f"{'--tcc ' if tcc else ''}"
            f"{f'--filter={filter} ' if filter else ''}"
            f"{f'--filter-threshold={filter_threshold} ' if filter_threshold else ''}"
            f"{'--overwrite ' if overwrite else ''}"
            f"{'--dry-run ' if dry_run else ''}"
            f"{'--loom ' if loom else ''}"
            f"{'--h5ad ' if h5ad else ''}"
            f"{'--cellranger ' if cellranger else ''}"
            f"{'--gene-names ' if gene_names else ''}"
            f"{'--report ' if report else ''}"
            f"{f'--kallisto={kallisto} ' if kallisto else ''}"
            f"{f'--bustools={bustools} ' if bustools else ''}"
            f"{f'-c1 {c1} ' if c1 else ''}"
            f"{f'-c2 {c2} ' if c2 else ''}"
            f"{f'--parity={parity} ' if parity else ''}"
            f"{f'--fragment-l={fragment_l} ' if fragment_l else ''}"
            f"{f'--fragment-s={fragment_s} ' if fragment_s else ''}"
            f"fastqs {files}"
        )

        os.system(
            command
        )

    def list(self):
        os.system(
            "kb list"
        )


