from pykallisto import Kallisto, KallistoBus
import pathlib 
import os 

here = pathlib.Path(__file__).parent.resolve()

k = KallistoBus(
    files=os.path.join(here, '..', 'local', 'pbmc_1k_v3_fastqs')
)

# k.ref(index='index', t2g='t2g', )