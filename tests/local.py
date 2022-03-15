from pykallisto import Kallisto 
import pathlib 
import os 

here = pathlib.Path(__file__).parent.resolve()

k = Kallisto(
    files=os.path.join(here, '..', 'local', 'fastq')
)
print(k.files)

k.index(index='index_file', threads=4)