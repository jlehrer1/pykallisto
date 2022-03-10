from pykallisto import Kallisto 

test = Kallisto(files=['asdf'])
print(test.quant(output_dir='asdf').strip())