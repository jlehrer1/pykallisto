from multiprocessing.sharedctypes import Value
import os 
import pathlib 
from typing import *

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
        output_path: str, 
        files: List[str]=None,
        
    ):
        if not files and not self.files:
            raise ValueError("One of self.files or files argument must not be None.")

        if files and self.files:
            print('Warning: Class initialized with files list & files list passed. Using passed list.')

        if not files:
            files = self.files

        if not self._index_file:
            raise ValueError("Index files not generated. Run self.index().")

        os.system(
            f"kallisto quant --index={self._index_file} --output={output_path} {''.join(*files)}"
        )


        