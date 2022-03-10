from multiprocessing.sharedctypes import Value
import os 
import pathlib 
from typing import *
import subprocess 

class Kallisto:
    def __init__(
        self, 
        files: List[str]=None,
    ):
        self.files = files 
        self._index_file = None 

    def index(
        self,
        index,
        kmer_size: int=None,
    ):  
        """Builds pseudoindex for Kalliso"""
        if not kmer_size:
            kmer_size = 31

        self._index_file = index 

        os.system(
            f"kallisto index --index={index} --kmer-size={kmer_size}"
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
        filestr = ''.join(*files)

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

        return command 

    def bus(
        self,
        output_dir: str,
        files: List[str],
        technology: str,
        list: bool=False,
        threads: int=1
    ) -> None: 
    
        os.system(
            f"kallisto bus --index={self._index.file} "
            f"--output-dir={output_dir}"
            f"--technology={technology}"
            f"--threads={threads}"
            f"{'--list ' if list else ''}"
            f"{''.join(*files)}"
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

        os.system(
            "kallisto pseudo "
            f"--index={self._index_file} ",
            f"--output-dir={output_dir} ",
            f"{'--umi ' if umi else ''}",
            f"{f'--batch={batch} ' if batch else ''}"
            f"{'--single ' if single else ''}"
            f"{f'--fragment-length={fragment_length} ' if fragment_length else ''}"
            f"{f'--sd={sd} ' if sd else ''}"
            f"--threads={threads}"
        )

class KB:
    def __init__(
        self,
        files: List[str]
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
        fasta: str,
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
        bustools: str=None
    ) -> None:
        pass 

    def count(
        self,
        files: List[str],
        index: str,
        t2g: str,
        technology: str,
        temp,
        keep_tmp,
        verbose,
        out,
        whilelist,
        threads,
        memory,
        strand,
        workflow,
        em,
        umi_gene,
        mm,
        tcc,
        filter,
        filter_threshold,
        overwrite,
        dry_run,
        loom,
        h5ad,
        cellranger,
        gene_names,
        report,
        kallisto,
        bustools,
        c1,
        c2,
        parity,
        fragment_l,
        fragment_s,
    ) -> None:
        pass 

    def 

