from pykallisto import Kallisto 
import pathlib 
import os 

here = pathlib.Path(__file__).parent.resolve()

k = Kallisto(
    files=os.path.join(here, '..', 'local', 'pbmc_1k_v3_fastqs')
)

k.index(index='index_file', make_unique=True)