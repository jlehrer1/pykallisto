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
        if isinstance(files, str):
            self.files = [os.path.join(files, f) for f in os.listdir(files)]
        elif isinstance(files, list):
            self.files = files 
        else:
            raise ValueError("files must be path to files or list of files")
        
        self._index_file = (index_file if index_file else None) 

    def index(
        self,
        index,
        kmer_size: int=None,
        make_unique: bool=False 
    ):  
        """Builds pseudoindex for Kalliso"""
        self._index_file = index 
        files = ' '.join(self.files)

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
        files: List[str]=None,
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
        chromosomes: str=None
    ):
        if not files and not self.files:
            raise ValueError("One of self.files or files argument must not be None.")

        if files and self.files:
            print('Warning: Class initialized with files list & files list passed. Using passed list.')

        if not files:
            files = self.files

        if not self._index_file:
            raise ValueError("Index files not generated. Run self.index().")
        filestr = ' '.join(files)

        command = (
            f"kallisto quant " 
            f"--index={self._index_file} "
            f"{'--bias ' if bias else ''}"
            f"{'--bootstrap_samples={bootstrap_samples} ' if bootstrap_samples else ''}"
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
            f"--output={output_dir} {filestr}"
        )

        os.system(
            command 
        ) 

    def bus(
        self,
        output_dir: str,
        files: List[str],
        technology: str,
        list: bool=False,
        threads: int=1
    ) -> None:

        files = ' '.join(files)
        command = (
            f"kallisto bus --index={self._index.file} "
            f"--output-dir={output_dir}"
            f"--technology={technology}"
            f"--threads={threads}"
            f"{'--list ' if list else ''}"
            f"{files}"
        )

        os.system(
            command 
        )

    def pseudo(
        self, 
        output_dir: str,
        umi: bool=False,
        batch: str=None,
        single: bool=False,
        fragment_length: float=None,
        sd: float=None,
        threads: int=1,
    ) -> None:

        command = ("kallisto pseudo "
            f"--index={self._index_file} ",
            f"--output-dir={output_dir} ",
            f"{'--umi ' if umi else ''}",
            f"{f'--batch={batch} ' if batch else ''}"
            f"{'--single ' if single else ''}"
            f"{f'--fragment-length={fragment_length} ' if fragment_length else ''}"
            f"{f'--sd={sd} ' if sd else ''}"
            f"--threads={threads}"
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
        
        command = (
            "kb ref "
            f"-i={index} "
            f"-g={t2g} "
            f"{f'-f1={fasta} ' if fasta else ''}"
            f"{'--tmp ' if tmp else ''}"
            f"{'--keep-tmp ' if keep_temp else ''}"
            f"{'--verbose ' if verbose else ''}"
            f"{include_attribute_str if include_attribute else ''}"
            f"{exclude_attribute_str if exclude_attribute else ''}"
            f"{f'-d={d} ' if d else ''}"
            f"{f'k={k} ' if k else ''}"
            f"{f'workflow={workflow} ' if workflow else ''}"
            f"{'--overwrite ' if overwrite else ''}"
            f"{f'--kallisto={kallisto} ' if kallisto else ''}"
            f"{f'--bustools={bustools} ' if bustools else ''}"
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


