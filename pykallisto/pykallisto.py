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
        index_file: str,
        output_dir: str,
        threads: int=1,
    ) -> None:

        os.system(
            ""
        )

class KB:
    def __init__(
        self,
        files,
    ) -> None:
        self.files = files 

    def info(self):
        pass 

    def ref(self, index_file):
        pass 
          